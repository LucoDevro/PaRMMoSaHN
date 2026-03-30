import logging
import shutil
from pathlib import Path

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


def run_pan(
    genomes_dir: Path,
    output_dir: Path,
    label: str,
    threads: int = 14,
    add_unique: bool = False,
):
    """Executes Step 1 & 2: Pangenome and Pan-model."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- STEP 1: Pangenome ---
    logger.info("=== STEP 1: Building Pangenome ===")
    pangenome_out = output_dir / "01-pangenome"
    softcore_fasta = build_pangenome(
        genomes_dir, pangenome_out, label, threads=threads, add_unique=add_unique
    )

    # --- STEP 2: Pan-model ---
    logger.info("=== STEP 2: Building Pan-model ===")
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

    # --- STEP 4: Automatic Curation ---
    curated_out = output_dir / "04-curated_models"
    curated_out.mkdir(exist_ok=True)

    if curation_db and curation_db.exists():
        logger.info("=== STEP 4: Applying Automatic Curations ===")
        for draft_model in draft_models_dir.glob("*/*.xml"):
            apply_curations(draft_model, curation_db, curated_out)
    else:
        logger.info("=== STEP 4: Automatic Curation (SKIPPED - No DB provided) ===")
        # Safely copy draft models forward so the pipeline doesn't break
        for draft_model in draft_models_dir.glob("*/*.xml"):
            shutil.copy(draft_model, curated_out / draft_model.name)

    current_models_dir = curated_out

    # --- STEP 5: ModelPolisher via Docker (Optional) ---
    if no_polish:
        logger.info("=== STEP 5: ModelPolisher (SKIPPED by user) ===")
    elif not shutil.which(engine if engine != "podman" else "podman-compose"):
        logger.warning(
            f"=== STEP 5: ModelPolisher (SKIPPED - {engine} is not installed) ==="
        )
    else:
        logger.info(f"=== STEP 5: ModelPolisher ({engine}) ===")
        polished_out = output_dir / "05-polished"
        current_models_dir = run_model_polisher(
            current_models_dir, polished_out, engine=engine
        )

    # --- STEP 6: Gather & Convert ---
    logger.info("=== STEP 6: Gathering Final Results ===")
    final_out = output_dir / "06-final_results"
    gather_and_convert_models(current_models_dir, final_out)

    # --- STEP 7: Final Evaluation ---
    evaluate_final_models(current_models_dir, output_dir)

    logger.info(f"Pipeline complete! All results saved to {final_out}")
