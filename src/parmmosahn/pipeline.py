import json
import logging
import shutil
import time
from datetime import datetime, timezone
from importlib.metadata import version as pkg_version
from pathlib import Path

from parmmosahn.annotation import enrich_models_directory
from parmmosahn.constants import step_banner
from parmmosahn.curation import (
    apply_curations,
    evaluate_final_models,
    generate_curation_template,
    run_model_polisher,
)
from parmmosahn.models import build_pan_model, derive_strain_models
from parmmosahn.pangenome import build_pangenome
from parmmosahn.results import gather_and_convert_models

logger = logging.getLogger(__name__)


def _format_elapsed(seconds: float) -> str:
    """Format elapsed seconds into a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s}s"
    h, m = divmod(m, 60)
    return f"{h}h {m}m {s}s"


def _count_xml(directory: Path) -> int:
    """Count XML files recursively in a directory."""
    if not directory.exists():
        return 0
    return len(list(directory.rglob("*.xml")))


def _collect_model_statistics(output_dir: Path) -> dict:
    """Scan final models and collect biological statistics for the summary."""
    import cobra as cb

    stats: dict = {}
    final_dir = output_dir / "06-final_models"
    if not final_dir.exists():
        return stats

    xml_files = sorted(final_dir.rglob("*.xml"))
    if not xml_files:
        return stats

    # Suppress noisy cobra INFO messages (e.g. "Compartment `e0` sounds like
    # an external compartment") during batch model loading.
    cobra_logger = logging.getLogger("cobra")
    prev_level = cobra_logger.level
    cobra_logger.setLevel(logging.WARNING)

    per_model: list[dict] = []
    all_rxns: set[str] = set()
    rxn_sets: list[set[str]] = []

    for xml in xml_files:
        try:
            model = cb.io.read_sbml_model(str(xml))
            rxns = {r.id for r in model.reactions}
            all_rxns.update(rxns)
            rxn_sets.append(rxns)
            per_model.append(
                {
                    "strain": xml.stem,
                    "reactions": len(model.reactions),
                    "metabolites": len(model.metabolites),
                    "genes": len(model.genes),
                    "exchanges": len(model.exchanges),
                }
            )
        except Exception:
            continue

    cobra_logger.setLevel(prev_level)

    if not per_model:
        return stats

    n = len(per_model)
    for key in ("reactions", "metabolites", "genes"):
        vals = [m[key] for m in per_model]
        mean = sum(vals) / n
        std = (sum((v - mean) ** 2 for v in vals) / n) ** 0.5
        stats[key] = {
            "min": min(vals),
            "max": max(vals),
            "mean": round(mean, 1),
            "std": round(std, 1),
        }

    # Core/accessory/unique reaction counts
    core = sum(1 for r in all_rxns if all(r in s for s in rxn_sets))
    unique = sum(1 for r in all_rxns if sum(1 for s in rxn_sets if r in s) == 1)
    stats["pan_reactome"] = {
        "total": len(all_rxns),
        "core": core,
        "accessory": len(all_rxns) - core - unique,
        "unique": unique,
    }

    stats["per_model"] = per_model
    return stats


def write_pipeline_summary(
    output_dir: Path,
    steps: dict,
    parameters: dict,
    command: str = "",
    started_at: str = "",
    total_elapsed: float = 0.0,
) -> Path:
    """Write a pipeline_summary.json file summarizing the run."""
    summary = {
        "parmmosahn_version": pkg_version("parmmosahn"),
        "command": command,
        "started_at": started_at,
        "finished_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "elapsed": _format_elapsed(total_elapsed),
        "steps": steps,
        "parameters": parameters,
    }

    # Collect biological statistics from final models
    model_stats = _collect_model_statistics(output_dir)
    if model_stats:
        summary["model_statistics"] = model_stats

    summary_path = output_dir / "pipeline_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info(f"Pipeline summary written to {summary_path}")
    return summary_path


def run_pan(
    genomes_dir: Path,
    output_dir: Path,
    label: str,
    threads: int = 14,
    add_unique: bool = False,
    soft_core_override: float | None = None,
):
    """Executes Step 1 & 2: Pangenome and Pan-model."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- STEP 1: Pangenome ---
    logger.info(step_banner(1))
    pangenome_out = output_dir / "01-pangenome"
    softcore_fasta = build_pangenome(
        genomes_dir,
        pangenome_out,
        label,
        threads=threads,
        add_unique=add_unique,
        soft_core_override=soft_core_override,
    )

    # --- STEP 2: Pan-model ---
    logger.info(step_banner(2))
    panmodel_out = output_dir / "02-panmodel"
    build_pan_model(softcore_fasta, panmodel_out, threads=threads)

    logger.info(f"Pan-tasks complete! Pan-model saved to {panmodel_out}")
    return panmodel_out / f"{label}-all-Reactions.tbl", softcore_fasta


def run_downstream(
    draft_models_dir: Path,
    output_dir: Path,
    curation_db: Path | None = None,
    no_polish: bool = False,
    engine: str = "docker",
):
    """Executes Steps 4-7: Curation, Polishing, and Export."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- STEP 4a: Generate curation template if none provided ---
    if curation_db and Path(curation_db).exists():
        logger.info("Using provided curation database: %s", curation_db)
    else:
        logger.info("=== STEP 4a: Generating curation template (memote on drafts) ===")
        curation_db = generate_curation_template(draft_models_dir, output_dir)
        logger.info(f"Curation template generated: {curation_db}")

    # --- STEP 4b: Automatic Curation ---
    curated_out = output_dir / "04-curated_models"
    curated_out.mkdir(exist_ok=True)

    logger.info(step_banner(4))
    for draft_model in draft_models_dir.rglob("*.xml"):
        apply_curations(draft_model, curation_db, curated_out)

    current_models_dir = curated_out

    # --- STEP 5: ModelPolisher via Docker (Optional) ---
    # Map engine name to the actual binary used by run_model_polisher so the
    # availability check reflects the real command (e.g. podman -> podman-compose).
    engine_binary_map = {
        "docker": "docker",
        "podman": "podman-compose",
        "apptainer": "apptainer",
    }
    engine_bin = engine_binary_map.get(engine, engine)
    if no_polish:
        logger.info(step_banner(5, "SKIPPED by user"))
    elif not shutil.which(engine_bin):
        logger.warning(step_banner(5, f"SKIPPED - {engine_bin} is not installed"))
    else:
        logger.info(step_banner(5, engine))
        polished_out = output_dir / "05-polished_models"
        polished_result = run_model_polisher(
            current_models_dir, polished_out, engine=engine
        )
        if _count_xml(polished_result) > 0:
            current_models_dir = polished_result
        else:
            logger.warning(
                "ModelPolisher produced no output — continuing with pre-polish models."
            )

    # --- STEP 5b: Annotation Enrichment ---
    logger.info("=== STEP 5b: Annotation Enrichment ===")
    annotated_out = output_dir / "05b-annotated_models"
    # Try to locate gene ID mapping files for Bakta→NCBI resolution
    proteomes_dir = draft_models_dir.parent / "proteomes"
    current_models_dir = enrich_models_directory(
        current_models_dir,
        annotated_out,
        proteomes_dir=proteomes_dir if proteomes_dir.is_dir() else None,
    )

    # --- STEP 6: Gather & Convert ---
    logger.info(step_banner(6))
    final_out = output_dir / "06-final_models"
    gather_and_convert_models(current_models_dir, final_out)

    # --- STEP 7: Final Evaluation ---
    evaluate_final_models(current_models_dir, output_dir)

    logger.info(f"Pipeline complete! All results saved to {final_out}")
