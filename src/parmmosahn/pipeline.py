import json
import logging
import shutil
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
    write_curation_application_report,
)
from parmmosahn.models import build_pan_model
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
    n_failed = 0

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
        except Exception as exc:
            n_failed += 1
            logger.warning(f"Statistics: failed to load {xml.name}: {exc}")

    cobra_logger.setLevel(prev_level)

    if n_failed:
        logger.warning(
            f"Statistics computed over {len(per_model)}/{len(xml_files)} final "
            f"models ({n_failed} failed to load)."
        )

    if not per_model:
        return stats

    # Record how many of the expected models actually contributed, so the
    # summary's biological statistics are never silently computed on a subset.
    stats["models_loaded"] = len(per_model)
    stats["models_expected"] = len(xml_files)

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


def _collect_provenance(output_dir: Path, parameters: dict) -> dict:
    """Build a per-run provenance block from parameters and on-disk artefacts.

    Captures:
      * SHA-256 of the medium CSV and curation database (so re-runs can be
        verified bit-for-bit against the original inputs).
      * The occurrence threshold actually used by PPanGGOLiN, recovered from
        the ``01-pangenome/pangenome_meta.json`` sidecar that
        :func:`build_pangenome` writes.

    Every field is best-effort — missing inputs are simply omitted so the
    summary survives runs that skip stages (e.g. ``pan`` only).
    """
    from parmmosahn.utils import sha256_file

    prov: dict = {}

    medium = parameters.get("medium")
    if medium:
        medium_path = Path(medium)
        if medium_path.exists():
            prov["medium_path"] = str(medium_path)
            prov["medium_sha256"] = sha256_file(medium_path)

    curation_db = parameters.get("curation_db")
    if curation_db:
        db_path = Path(curation_db)
        if db_path.exists():
            prov["curation_db_path"] = str(db_path)
            prov["curation_db_sha256"] = sha256_file(db_path)

    meta_path = output_dir / "01-pangenome" / "pangenome_meta.json"
    if meta_path.exists():
        try:
            pan_meta = json.loads(meta_path.read_text(encoding="utf-8"))
            prov["pangenome"] = pan_meta
        except (json.JSONDecodeError, OSError):
            pass

    return prov


def write_pipeline_summary(
    output_dir: Path,
    steps: dict,
    parameters: dict,
    command: str = "",
    started_at: str = "",
    total_elapsed: float = 0.0,
) -> Path:
    """Write a pipeline_summary.json file summarizing the run.

    Merges with any existing summary so that resumed or re-run pipelines
    preserve prior step timings and parameters instead of overwriting them.
    """
    summary_path = output_dir / "pipeline_summary.json"

    # Load existing summary to preserve prior run data
    existing: dict = {}
    if summary_path.exists():
        try:
            existing = json.loads(summary_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            pass

    # Merge steps: prior steps preserved, new/updated steps take precedence
    merged_steps = existing.get("steps", {})
    merged_steps.update(steps)

    # Merge parameters: prior values preserved, new ones added/updated
    merged_params = existing.get("parameters", {})
    merged_params.update(parameters)

    # Capture external tool versions for reproducibility. Probe once per
    # write and merge with any previously recorded versions so the entry
    # survives across resumed runs.
    from parmmosahn.utils import check_environment, get_system_resources

    env_probe = check_environment()
    external_versions = existing.get("external_tool_versions", {})
    for tool, info in env_probe.items():
        if info.get("found") and info.get("version"):
            external_versions[tool] = info["version"]

    # Run environment is captured once and pinned to the first successful
    # summary write; later writes only add fields, not overwrite them, so
    # a resumed run preserves the host where the bulk of work was done.
    run_environment = existing.get("run_environment") or get_system_resources()

    provenance = existing.get("provenance", {})
    provenance.update(_collect_provenance(output_dir, merged_params))

    summary = {
        "parmmosahn_version": pkg_version("parmmosahn"),
        "external_tool_versions": external_versions,
        "command": command,
        "started_at": existing.get("started_at") or started_at,
        "finished_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "elapsed": _format_elapsed(total_elapsed),
        "run_environment": run_environment,
        "provenance": provenance,
        "steps": merged_steps,
        "parameters": merged_params,
    }

    # Collect biological statistics from final models
    model_stats = _collect_model_statistics(output_dir)
    if model_stats:
        summary["model_statistics"] = model_stats

    # Atomic write: this file is read back on the next run to merge prior
    # provenance/timings, so a truncated write would silently wipe run history.
    from parmmosahn.utils import atomic_write_text

    atomic_write_text(summary_path, json.dumps(summary, indent=2))
    logger.info(f"Pipeline summary written to {summary_path}")
    return summary_path


def run_pan(
    genomes_dir: Path,
    output_dir: Path,
    label: str,
    threads: int = 14,
    add_small: bool = False,
    occurrence_override: float | None = None,
):
    """Executes Step 1 & 2: Pangenome and Pan-model."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- STEP 1: Pangenome ---
    logger.info(step_banner(1))
    pangenome_out = output_dir / "01-pangenome"
    included_fasta = build_pangenome(
        genomes_dir,
        pangenome_out,
        label,
        threads=threads,
        add_small=add_small,
        occurrence_override=occurrence_override,
    )

    # --- STEP 2: Pan-model ---
    logger.info(step_banner(2))
    panmodel_out = output_dir / "02-panmodel"
    build_pan_model(included_fasta, panmodel_out, threads=threads)

    logger.info(f"Pan-tasks complete! Pan-model saved to {panmodel_out}")
    return panmodel_out / f"{label}-all-Reactions.tbl", included_fasta
