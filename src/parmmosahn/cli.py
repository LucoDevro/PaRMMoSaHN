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


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option()
def cli():
    """
    PaRMMoSaHN: Pangenome Reference-based Metabolic Modelling by Saving Homolog Networks.
    """
    pass


@cli.command()
@click.option(
    "-g",
    "--genomes",
    required=True,
    type=click.Path(exists=True),
    help="Directory of input genome files (.gbff format).",
)
@click.option(
    "-o",
    "--output",
    default="01-pangenome",
    type=click.Path(),
    show_default=True,
    help="Target directory for all pangenome output files.",
)
@click.option(
    "-l",
    "--label",
    required=True,
    help="A descriptive label for the pangenome (e.g., 'clostridiales').",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
def pangenome(genomes, output, label, threads):
    """Step 1: Construct the partitioned pangenome using PPanGGOLiN."""
    click.echo(f"Building pangenome '{label}'...")
    build_pangenome(genomes, output, label, threads=threads)


@cli.command()
@click.option(
    "-f",
    "--fasta",
    required=True,
    type=click.Path(exists=True),
    help="Path to the PPanGGOLiN soft-core FASTA file.",
)
@click.option(
    "-o",
    "--output",
    default="02-pan_model",
    type=click.Path(),
    show_default=True,
    help="Target directory for the pan-model output.",
)
@click.option(
    "-t",
    "--threads",
    default=14,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
def panmodel(fasta, output, threads):
    """Step 2: Build the species-wide pan-model using gapseq."""
    click.echo("Building pan-model...")
    build_pan_model(fasta, output, threads=threads)


@cli.command()
@click.option(
    "-g",
    "--genomes",
    required=True,
    type=click.Path(exists=True),
    help="Directory of strain genome files (.gbff or .faa format).",
)
@click.option(
    "-f",
    "--fasta",
    required=True,
    type=click.Path(exists=True),
    help="Path to the PPanGGOLiN soft-core FASTA file.",
)
@click.option(
    "-r",
    "--panmodel-rxns",
    required=True,
    type=click.Path(exists=True),
    help="Path to the pan-model reaction table file.",
)
@click.option(
    "-o",
    "--output",
    default="03-strain_models",
    type=click.Path(),
    show_default=True,
    help="Target directory for the generated strain-specific draft models.",
)
@click.option(
    "-m",
    "--medium",
    required=True,
    type=click.Path(exists=True),
    help="Path to the media CSV file used for gapfilling.",
)
@click.option(
    "-t",
    "--threads",
    default=14,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
def strains(genomes, fasta, panmodel_rxns, output, medium, threads):
    """Step 3: Derive strain-specific models using diamond and gapseq."""
    click.echo("Deriving strain models...")
    derive_strain_models(genomes, fasta, panmodel_rxns, output, medium, threads=threads)


@cli.command()
@click.option(
    "--source-model",
    required=True,
    type=click.Path(exists=True),
    help="Path to the input draft SBML model (.xml).",
)
@click.option(
    "-c",
    "--curated-model",
    required=True,
    type=click.Path(),
    help="Output path for the curated SBML model.",
)
@click.option(
    "--db",
    required=True,
    type=click.Path(exists=True),
    help="Path to the Excel Curation Database containing curated reactions.",
)
def curate(source_model, curated_model, db):
    """Step 4: Apply curation database rules to a specific model."""
    click.echo(f"Curating model {Path(source_model).name}...")
    apply_curations(Path(source_model), Path(curated_model), Path(db))


@cli.command()
@click.option(
    "-M",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing the curated SBML models.",
)
@click.option(
    "-o",
    "--output",
    default="05-polished_models",
    type=click.Path(),
    show_default=True,
    help="Target directory for the polished models.",
)
@click.option(
    "-e",
    "--engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
    show_default=True,
    help="Container engine to employ for running ModelPolisher.",
)
def polish(models, output, engine):
    """Step 5: Polish models using Palsson lab ModelPolisher (requires a container engine)."""
    from parmmosahn.curation import run_model_polisher

    run_model_polisher(models, output, engine=engine)


@cli.command()
@click.option(
    "-M",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing the final SBML models.",
)
@click.option(
    "-o",
    "--output",
    default="06-final_models",
    type=click.Path(),
    show_default=True,
    help="Target directory for the exported output formats (JSON/MATLAB).",
)
def gather(models, output):
    """Step 6: Gather models and convert to JSON/MATLAB formats."""
    click.echo("Gathering and converting models...")
    gather_and_convert_models(models, output)


@cli.command()
@click.option(
    "-g",
    "--genomes",
    required=True,
    type=click.Path(exists=True),
    help="Directory of new strain genome files (.faa or .gbff format).",
)
@click.option(
    "-f",
    "--fasta",
    required=True,
    type=click.Path(exists=True),
    help="Path to an existing soft-core FASTA file.",
)
@click.option(
    "-r",
    "--panmodel-rxns",
    required=True,
    type=click.Path(exists=True),
    help="Path to an existing pan-model reaction table file.",
)
@click.option(
    "-m",
    "--medium",
    required=True,
    type=click.Path(exists=True),
    help="Path to the media CSV file for gapfilling.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Target output directory for the new models.",
)
@click.option(
    "-c",
    "--curation-db",
    type=click.Path(exists=True),
    help="Optional path to an Excel DB for automated curation.",
)
@click.option(
    "-t",
    "--threads",
    default=14,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    "--skip-polisher", is_flag=True, help="Skip the ModelPolisher execution step."
)
@click.option(
    "-e",
    "--engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
    show_default=True,
    help="Container engine to use for running tools like ModelPolisher.",
)
@click.option(
    "--stop-after-drafts",
    is_flag=True,
    help="Pause the pipeline before automatic curation to allow manual evaluation/edits.",
)
def project(
    genomes,
    fasta,
    panmodel_rxns,
    medium,
    output,
    curation_db,
    threads,
    skip_polisher,
    engine,
    stop_after_drafts,
):
    """
    Rapidly project new strains (.faa or .gbff format) onto an existing pan-model and continue with downstream curation.
    This bypasses Steps 1 & 2 (PPanGGOLiN).
    """
    from parmmosahn.curation import generate_curation_template
    from parmmosahn.models import derive_strain_models
    from parmmosahn.pipeline import run_downstream

    out_path = Path(output)
    out_path.mkdir(parents=True, exist_ok=True)

    # 1. Run Strain Derivation (Step 3)
    click.echo("=== STEP 3: Deriving Strain Models ===")
    strains_out = out_path / "03-strains"
    draft_models_dir = derive_strain_models(
        genomes_dir=genomes,
        pangenome_fasta=fasta,
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
            f"When finished, resume using: parmmosahn resume -o {out_path.resolve()}"
        )
        return

    # 2. Pipe directly into downstream processing (Steps 4-7)
    run_downstream(
        draft_models_dir=draft_models_dir,
        output_dir=out_path,
        curation_db=Path(curation_db) if curation_db else None,
        skip_polisher=skip_polisher,
        engine=engine,
    )


@cli.command()
@click.option(
    "-g",
    "--genomes",
    required=True,
    type=click.Path(exists=True),
    help="Directory of input genomes (.gbff format is required for the full pipeline, .faa format accepted for strain derivation).",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Target base directory for all outputs.",
)
@click.option(
    "-l", "--label", required=True, help="A descriptive label for the pangenome run."
)
@click.option(
    "-m",
    "--medium",
    required=True,
    type=click.Path(exists=True),
    help="Path to the media CSV file for gapfilling.",
)
@click.option(
    "-c",
    "--curation-db",
    required=False,
    type=click.Path(exists=True),
    help="Optional path to an Excel DB used for automated curation.",
)
@click.option(
    "-t",
    "--threads",
    default=14,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    "--skip-polisher", is_flag=True, help="Skip the ModelPolisher execution step."
)
@click.option(
    "--stop-after-drafts",
    is_flag=True,
    help="Pause the pipeline right before curation to allow manual evaluation/edits.",
)
@click.option(
    "-e",
    "--engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
    show_default=True,
    help="Container engine to use for running ModelPolisher.",
)
def run(
    genomes,
    output,
    label,
    medium,
    curation_db,
    threads,
    skip_polisher,
    stop_after_drafts,
    engine,
):
    """Run the complete PaRMMoSaHN pipeline from start to finish."""
    from parmmosahn.pipeline import run_all

    run_all(
        genomes_dir=Path(genomes),
        output_dir=Path(output),
        label=label,
        medium_file=Path(medium),
        curation_db=Path(curation_db) if curation_db else None,
        threads=threads,
        skip_polisher=skip_polisher,
        engine=engine,
        stop_after_drafts=stop_after_drafts,
    )


@cli.command()
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Base output directory from the previously paused run.",
)
@click.option(
    "-d",
    "--drafts",
    type=click.Path(exists=True, path_type=Path),
    help="Optional path: Folder containing draft models. Defaults to <output>/03-strains/models.",
)
@click.option(
    "-c",
    "--curation-db",
    type=click.Path(path_type=Path),
    help="Optional path: Path to the Excel DB. Defaults to <output>/curation_template.xlsx.",
)
@click.option(
    "--skip-polisher", is_flag=True, help="Skip the ModelPolisher execution step."
)
@click.option(
    "-e",
    "--engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
    show_default=True,
    help="Container engine to use for running ModelPolisher.",
)
def resume(
    output: Path,
    drafts: Path | None,
    curation_db: Path | None,
    skip_polisher,
    engine,
):
    """Resume the configured pipeline originating from Step 4 using a finished curation database."""
    from parmmosahn.pipeline import run_downstream

    # Auto-resolve defaults based on the project root
    if not drafts:
        drafts = output / "03-strains" / "models"
        if not drafts.exists():
            raise click.BadParameter(
                f"Could not find draft models at {drafts}. Please specify --drafts manually."
            )

    if not curation_db:
        curation_db = output / "curation_template.xlsx"

    # Informative error if the default file isn't there
    if not curation_db.exists():
        raise click.FileError(
            str(curation_db),
            hint="Could not find 'curation_template.xlsx' in the designated output directory. "
            "If you renamed or moved your compiled Excel file, please define its exact path via the --curation-db flag.",
        )

    click.echo(f"Resuming pipeline referencing curation database: {curation_db.name}")

    run_downstream(
        draft_models_dir=drafts,
        output_dir=output,
        curation_db=curation_db,
        skip_polisher=skip_polisher,
        engine=engine,
    )


def main():
    cli()


if __name__ == "__main__":
    main()
