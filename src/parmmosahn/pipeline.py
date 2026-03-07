import logging
import shutil
from pathlib import Path

from parmmosahn.curation import apply_curations, run_model_polisher, run_sbo_annotation
from parmmosahn.models import build_pan_model, derive_strain_models
from parmmosahn.pangenome import build_pangenome
from parmmosahn.results import gather_and_convert_models

logger = logging.getLogger(__name__)


def run_all(
    genomes_dir: Path,
    output_dir: Path,
    label: str,
    medium_file: Path,
    curation_db: Path,
    threads: int = 14,
    sbo_wrapper: Path | None = None,
    sbo_dbs: Path | None = None,
    skip_polisher: bool = False,
    engine: str = "docker",
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

    # --- STEP 4: Automatic Curation ---
    logger.info("=== STEP 4: Applying Automatic Curations ===")
    curated_out = output_dir / "04-curated_models"
    curated_out.mkdir(exist_ok=True)

    for draft_model in draft_models_dir.glob("*/*.xml"):
        apply_curations(draft_model, curation_db, curated_out)

    current_models_dir = curated_out

    # --- STEP 5: SBO Annotation (Optional) ---
    if sbo_wrapper and sbo_dbs:
        logger.info("=== STEP 5: SBO Annotation ===")
        sbo_out = output_dir / "05-sbo_annotated"
        current_models_dir = run_sbo_annotation(
            current_models_dir, sbo_out, sbo_wrapper, sbo_dbs
        )
    else:
        logger.info("=== STEP 5: SBO Annotation (SKIPPED - Missing paths) ===")

    # --- STEP 6: ModelPolisher via Docker (Optional) ---
    if skip_polisher:
        logger.info("=== STEP 6: ModelPolisher (SKIPPED by user) ===")
    elif not shutil.which(
        engine if engine != "podman" else "podman-compose"
    ):  # Quick check if tool exists
        logger.warning(
            f"=== STEP 6: ModelPolisher (SKIPPED - {engine} is not installed) ==="
        )
    else:
        logger.info(f"=== STEP 6: ModelPolisher ({engine}) ===")
        polished_out = output_dir / "06-polished"
        current_models_dir = run_model_polisher(
            current_models_dir, polished_out, engine=engine
        )

    # --- STEP 7: Gather & Convert ---
    logger.info("=== STEP 7: Gathering Final Results ===")
    final_out = output_dir / "07-final_results"
    gather_and_convert_models(current_models_dir, final_out)

    logger.info(f"Pipeline complete! All results saved to {final_out}")
