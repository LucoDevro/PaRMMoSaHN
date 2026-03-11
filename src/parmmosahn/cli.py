from pathlib import Path

import click

from parmmosahn.curation import (
    apply_curations,
    extract_memote_problems,
    run_memote_snapshot,
)
from parmmosahn.models import build_pan_model, derive_strain_models

# Import our functional modules
from parmmosahn.pangenome import build_pangenome
from parmmosahn.results import gather_and_convert_models


@click.group()
@click.version_option()
def cli():
    """
    PaRMMoSaHN: Pangenome Reference-based Metabolic Modelling by Saving Homolog Networks.
    """
    pass


@cli.command()
@click.option(
    "--genomes-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory of input .gbff genomes.",
)
@click.option(
    "--output-dir", default="01-pangenome", type=click.Path(), help="Output directory."
)
@click.option(
    "--label", required=True, help="Label for the pangenome (e.g., 'clostridiales')."
)
@click.option("--threads", default=1, help="Number of CPU threads.")
def pangenome(genomes_dir, output_dir, label, threads):
    """Step 1: Construct the partitioned pangenome using PPanGGOLiN."""
    click.echo(f"Building pangenome '{label}'...")
    build_pangenome(genomes_dir, output_dir, label, threads=threads)


@cli.command()
@click.option(
    "--fasta",
    required=True,
    type=click.Path(exists=True),
    help="Path to PPanGGOLiN soft-core FASTA.",
)
@click.option(
    "--output-dir", default="02-pan_model", type=click.Path(), help="Output directory."
)
@click.option("--threads", default=14, help="Number of CPU threads.")
def panmodel(fasta, output_dir, threads):
    """Step 2: Build the species-wide pan-model using gapseq."""
    click.echo("Building pan-model...")
    build_pan_model(fasta, output_dir, threads=threads)


@cli.command()
@click.option(
    "--genomes-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory of strain genomes (.gbff or .faa).",
)
@click.option(
    "--pangenome-fasta",
    required=True,
    type=click.Path(exists=True),
    help="Soft-core FASTA.",
)
@click.option(
    "--panmodel-rxns",
    required=True,
    type=click.Path(exists=True),
    help="Panmodel reaction table.",
)
@click.option(
    "--output-dir",
    default="03-strain_models",
    type=click.Path(),
    help="Output directory.",
)
@click.option(
    "--medium",
    required=True,
    type=click.Path(exists=True),
    help="Media CSV file for gapfilling.",
)
@click.option("--threads", default=14, help="Number of CPU threads.")
def strains(genomes_dir, pangenome_fasta, panmodel_rxns, output_dir, medium, threads):
    """Step 3: Derive strain-specific models using diamond and gapseq."""
    click.echo("Deriving strain models...")
    derive_strain_models(
        genomes_dir, pangenome_fasta, panmodel_rxns, output_dir, medium, threads=threads
    )


@cli.command()
@click.option(
    "--source-model",
    required=True,
    type=click.Path(exists=True),
    help="Input SBML model.",
)
@click.option(
    "--curated-model",
    required=True,
    type=click.Path(),
    help="Output curated SBML model.",
)
@click.option(
    "--db", required=True, type=click.Path(exists=True), help="Excel Curation Database."
)
def curate(source_model, curated_model, db):
    """Step 5: Apply curation database rules to a model."""
    click.echo(f"Curating model {Path(source_model).name}...")
    apply_curations(Path(source_model), Path(curated_model), Path(db))


@cli.command()
@click.option(
    "--model-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory of curated SBML models.",
)
@click.option(
    "--output-dir",
    default="06-final_models",
    type=click.Path(),
    help="Output directory.",
)
def gather(model_dir, output_dir):
    """Step 6: Gather models and convert to JSON/MATLAB."""
    click.echo("Gathering and converting models...")
    gather_and_convert_models(model_dir, output_dir)


@cli.command()
@click.option("--model-dir", required=True, type=click.Path(exists=True))
@click.option("--output-dir", default="04-annotated_models")
@click.option(
    "--sbo-script",
    required=True,
    type=click.Path(exists=True),
    help="Path to SBOannotator_wrapper.py",
)
@click.option(
    "--sbo-dbs", required=True, type=click.Path(), help="Path to SBO DBs directory"
)
def annotate_sbo(model_dir, output_dir, sbo_script, sbo_dbs):
    """Add SBO terms to curated models."""
    from parmmosahn.curation import run_sbo_annotation

    run_sbo_annotation(model_dir, output_dir, sbo_script, sbo_dbs)


@cli.command()
@click.option("--model-dir", required=True, type=click.Path(exists=True))
@click.option("--output-dir", default="05-polished_models")
def polish(model_dir, output_dir, polisher_dir):
    """Polish models using Palsson lab ModelPolisher (requires Docker)."""
    from parmmosahn.curation import run_model_polisher

    run_model_polisher(model_dir, output_dir, polisher_dir)


@cli.command()
@click.option(
    "--genomes-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory of new strain genomes (.faa or .gbff).",
)
@click.option(
    "--pangenome-fasta",
    required=True,
    type=click.Path(exists=True),
    help="Existing soft-core FASTA.",
)
@click.option(
    "--panmodel-rxns",
    required=True,
    type=click.Path(exists=True),
    help="Existing panmodel reaction table.",
)
@click.option(
    "--medium",
    required=True,
    type=click.Path(exists=True),
    help="Media CSV file for gapfilling.",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Output directory for new models.",
)
@click.option(
    "--curation-db",
    type=click.Path(exists=True),
    help="Optional Excel DB for curation.",
)
@click.option("--threads", default=14)
@click.option("--skip-polisher", is_flag=True, help="Skip the ModelPolisher step.")
@click.option(
    "--container-engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
)
@click.option(
    "--stop-after-drafts",
    is_flag=True,
    help="Stop pipeline before curation to allow manual evaluation.",
)
def project(
    genomes_dir,
    pangenome_fasta,
    panmodel_rxns,
    medium,
    output_dir,
    curation_db,
    threads,
    skip_polisher,
    container_engine,
    stop_after_drafts,
):
    """
    Rapidly project new strains (.faa or .gbff) onto an existing pan-model and run downstream curation.
    Skips Steps 1 & 2 (PPanGGOLiN).
    """
    from parmmosahn.curation import generate_curation_template
    from parmmosahn.models import derive_strain_models
    from parmmosahn.pipeline import run_downstream

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # 1. Run Strain Derivation (Step 3)
    click.echo("=== STEP 3: Deriving Strain Models ===")
    strains_out = out_path / "03-strains"
    draft_models_dir = derive_strain_models(
        genomes_dir=genomes_dir,
        pangenome_fasta=pangenome_fasta,
        panmodel_rxn_table=panmodel_rxns,
        output_dir=strains_out,
        medium_file=medium,
        threads=threads,
    )

    if stop_after_drafts:
        click.echo("=== Preparing Human-in-the-Loop Curation (running memote) ===")
        generate_curation_template(draft_models_dir, out_path)
        click.echo("=== Pipeline paused for Human-in-the-Loop Curation ===")
        click.echo(
            f"ACTION REQUIRED: Please fill out the curation template at {out_path}/curation_template.xlsx"
        )
        click.echo(
            f"When finished, resume using: parmmosahn resume --output-dir {out_path.resolve()}"
        )
        return

    # 2. Pipe directly into downstream processing (Steps 4-7)
    run_downstream(
        draft_models_dir=draft_models_dir,
        output_dir=out_path,
        curation_db=Path(curation_db) if curation_db else None,
        skip_polisher=skip_polisher,
        engine=container_engine,
    )


@cli.command()
@click.option(
    "--genomes-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory of input genomes (.gbff required for full pipeline, .faa accepted for strain derivation).",
)
@click.option("--output-dir", required=True, type=click.Path())
@click.option("--label", required=True)
@click.option("--medium", required=True, type=click.Path(exists=True))
@click.option(
    "--curation-db",
    required=False,
    type=click.Path(exists=True),
    help="Optional Excel DB for curation.",
)
@click.option("--threads", default=14)
@click.option("--skip-polisher", is_flag=True, help="Skip the ModelPolisher step.")
@click.option(
    "--stop-after-drafts",
    is_flag=True,
    help="Stop pipeline before curation to allow manual evaluation.",
)
@click.option(
    "--container-engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
    help="Container engine to use for ModelPolisher.",
)
def run(
    genomes_dir,
    output_dir,
    label,
    medium,
    curation_db,
    threads,
    skip_polisher,
    stop_after_drafts,
    container_engine,
):
    """Run the PaRMMoSaHN pipeline from start to finish."""
    from parmmosahn.pipeline import run_all

    run_all(
        genomes_dir=Path(genomes_dir),
        output_dir=Path(output_dir),
        label=label,
        medium_file=Path(medium),
        curation_db=Path(curation_db) if curation_db else None,
        threads=threads,
        skip_polisher=skip_polisher,
        engine=container_engine,
        stop_after_drafts=stop_after_drafts,
    )


@cli.command()
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Base output directory",
)
@click.option(
    "--drafts-dir",
    type=click.Path(exists=True, path_type=Path),
    help="Optional: Folder containing draft models. Defaults to <output-dir>/03-strains/models",
)
@click.option(
    "--curation-db",
    type=click.Path(path_type=Path),
    help="Optional: Path to Excel DB. Defaults to <output-dir>/curation_template.xlsx",
)
@click.option("--skip-polisher", is_flag=True, help="Skip the ModelPolisher step.")
@click.option(
    "--container-engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
)
def resume(
    output_dir: Path,
    drafts_dir: Path | None,
    curation_db: Path | None,
    skip_polisher,
    container_engine,
):
    """Resume the pipeline from Step 4 using a completed curation database."""
    from parmmosahn.pipeline import run_downstream

    # Auto-resolve defaults based on the project root
    if not drafts_dir:
        drafts_dir = output_dir / "03-strains" / "models"
        if not drafts_dir.exists():
            raise click.BadParameter(
                f"Could not find draft models at {drafts_dir}. Please specify --drafts-dir manually."
            )

    if not curation_db:
        curation_db = output_dir / "curation_template.xlsx"

    # Informative error if the default file isn't there
    if not curation_db.exists():
        raise click.FileError(
            str(curation_db),
            hint="Could not find 'curation_template.xlsx' in the output directory. "
            "If you moved or renamed your Excel file, please provide its exact path using the --curation-db flag.",
        )

    click.echo(f"Resuming pipeline using curation database: {curation_db.name}")

    run_downstream(
        draft_models_dir=drafts_dir,
        output_dir=output_dir,
        curation_db=curation_db,
        skip_polisher=skip_polisher,
        engine=container_engine,
    )


def main():
    cli()


if __name__ == "__main__":
    main()
