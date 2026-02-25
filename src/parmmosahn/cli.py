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
    help="Directory of strain .faa genomes.",
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
    apply_curations(source_model, curated_model, db)


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


def main():
    cli()


if __name__ == "__main__":
    main()
