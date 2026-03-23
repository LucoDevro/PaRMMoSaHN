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


def run_downstream(
    draft_models_dir: Path,
    output_dir: Path,
    curation_db: Path | None = None,
    skip_polisher: bool = False,
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
    if skip_polisher:
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
    # current_models_dir holds the final XMLs right before they were gathered
    evaluate_final_models(current_models_dir, output_dir)

    logger.info(f"Pipeline complete! All results saved to {final_out}")


def run_all(
    genomes_dir: Path,
    output_dir: Path,
    label: str,
    medium_file: Path,
    curation_db: Path | None = None,
    threads: int = 14,
    skip_polisher: bool = False,
    engine: str = "docker",
    stop_after_drafts: bool = False,
):
    """Executes the entire PaRMMoSaHN pipeline from end to end."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- STEP 1: Pangenome ---
    logger.info("=== STEP 1: Building Pangenome ===")
    pangenome_out = output_dir / "01-pangenome"
    softcore_fasta = build_pangenome(genomes_dir, pangenome_out, label, threads=threads)

    # --- STEP 2: Pan-model ---
    logger.info("=== STEP 2: Building Pan-model ===")
    panmodel_out = output_dir / "02-panmodel"
    build_pan_model(softcore_fasta, panmodel_out, threads=threads)
    panmodel_rxns = panmodel_out / f"{label}-all-Reactions.tbl"

    # --- STEP 3: Strain Derivation ---
    logger.info("=== STEP 3: Deriving Strain Models ===")
    strains_out = output_dir / "03-strains"
    draft_models_dir = derive_strain_models(
        genomes_dir,
        softcore_fasta,
        panmodel_rxns,
        strains_out,
        medium_file,
        threads=threads,
    )

    if stop_after_drafts:
        logger.info("=== Preparing Human-in-the-Loop Curation (running memote) ===")
        template_path = generate_curation_template(draft_models_dir, output_dir)
        logger.info("=== Pipeline paused for Human-in-the-Loop Curation ===")
        logger.info(f"Draft models are ready in: {draft_models_dir}")
        logger.info(
            f"ACTION REQUIRED: Please fill out the curation template at {template_path.resolve()}"
        )
        logger.info(
            f"When finished, resume using: parmmosahn resume -o {output_dir.resolve()}"
        )
        return

    # Pass everything to the downstream runner
    run_downstream(
        draft_models_dir=draft_models_dir,
        output_dir=output_dir,
        curation_db=curation_db,
        skip_polisher=skip_polisher,
        engine=engine,
    )
