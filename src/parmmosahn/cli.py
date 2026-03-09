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


@cli.group()
def evaluate():
    """Commands for evaluating draft models via Memote to generate curation TSVs."""
    pass


@evaluate.command()
@click.option(
    "--model",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to a single draft .xml model.",
)
@click.option(
    "--out-dir",
    required=True,
    type=click.Path(path_type=Path),
    help="Directory to save the .html and .json reports.",
)
def snapshot(model: Path, out_dir: Path):
    """Run memote on a single model to generate JSON and HTML reports."""
    from parmmosahn.curation import run_memote_snapshot

    out_json = out_dir / f"{model.stem}.json"
    out_html = out_dir / f"{model.stem}.html"

    click.echo(f"Running memote on {model.name}...")
    run_memote_snapshot(model, out_json, out_html)
    click.echo(f"Done! Reports saved to {out_dir}")


@evaluate.command()
@click.option(
    "--report-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Dir containing memote .json reports.",
)
@click.option(
    "--model-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Dir containing the drafted .xml models.",
)
@click.option(
    "--out-excel",
    default="curation_template.xlsx",
    show_default=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="File path to save the output Excel template.",
)
def extract(report_dir: Path, model_dir: Path, out_excel: Path):
    """Parse a directory of memote JSONs and extract imbalances/duplicates into an Excel template."""
    from parmmosahn.curation import extract_memote_problems

    click.echo(f"Extracting problems from {report_dir}...")

    extract_memote_problems(report_dir, model_dir, out_excel)

    # Using .resolve() prints the absolute path so the user knows exactly where it went!
    click.echo(f"Done! Excel template saved to {out_excel.resolve()}")


@cli.command()
@click.option("--genomes-dir", required=True, type=click.Path(exists=True))
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
    "--drafts-dir",
    required=True,
    type=click.Path(exists=True),
    help="Folder containing draft models (e.g., 03-strains)",
)
@click.option(
    "--output-dir", required=True, type=click.Path(), help="Base output directory"
)
@click.option(
    "--curation-db",
    required=True,
    type=click.Path(exists=True),
    help="Filled Excel DB for curation.",
)
@click.option("--skip-polisher", is_flag=True, help="Skip the ModelPolisher step.")
@click.option(
    "--container-engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
)
def resume(drafts_dir, output_dir, curation_db, skip_polisher, container_engine):
    """Resume the pipeline from Step 4 using a completed curation database."""
    from parmmosahn.pipeline import run_downstream

    run_downstream(
        draft_models_dir=Path(drafts_dir),
        output_dir=Path(output_dir),
        curation_db=Path(curation_db),
        skip_polisher=skip_polisher,
        engine=container_engine,
    )


def main():
    cli()


if __name__ == "__main__":
    main()
