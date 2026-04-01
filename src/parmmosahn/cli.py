import logging
from pathlib import Path

import click

from parmmosahn.constants import STEPS, step_banner
from parmmosahn.curation import apply_curations
from parmmosahn.models import build_pan_model, derive_strain_models
from parmmosahn.pangenome import build_pangenome
from parmmosahn.results import gather_and_convert_models

# Set up logging globally here in the CLI entrypoint
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


@click.group(
    context_settings=dict(
        help_option_names=["-h", "--help"],
        auto_envvar_prefix="PARMMOSAHN",
    )
)
@click.version_option()
def cli():
    """
    PaRMMoSaHN: Pangenome Reference-based Metabolic Modelling by Saving Homolog Networks.

    All options can also be set via environment variables prefixed with PARMMOSAHN_.
    For example: PARMMOSAHN_THREADS=32, PARMMOSAHN_BIOMASS=neg.
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
    default="00-pan_results",
    type=click.Path(),
    show_default=True,
    help="Target base directory for pan-related output files.",
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
    default=14,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    "--add-unique",
    is_flag=True,
    help="Whether to add unique genes (those not in soft-core) to the pan-model. This can increase coverage but may add noise.",
)
@click.option(
    "--soft-core",
    default=None,
    type=float,
    show_default=True,
    help=(
        "Override the soft-core frequency threshold (fraction 0–1). "
        "Defaults to 2/n_genomes, the biologically motivated formula. "
        "Only change this if you have a specific reason."
    ),
)
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True),
    help="Path to a YAML config file. CLI arguments override config values. Requires pyyaml (pip install 'parmmosahn[config]').",
)
def pan(genomes, output, label, threads, add_unique, soft_core, config):
    """Step 1 & 2: Construct the partitioned pangenome and build the pan-model."""
    from parmmosahn.pipeline import run_pan

    if config:
        from parmmosahn.utils import load_config

        cfg = load_config(Path(config))
        genomes = genomes or cfg.get("genomes")
        output = output or cfg.get("output", "00-pan_results")
        label = label or cfg.get("label")
        threads = threads if threads != 14 else cfg.get("threads", 14)
        add_unique = add_unique or cfg.get("add_unique", False)
        if soft_core is None:
            soft_core = cfg.get("soft_core")
    # Ensure the required values are available (may be provided via config)
    assert (
        genomes is not None and label is not None
    ), "Missing required arguments: genomes and label must be provided either via CLI or config"

    click.echo(f"Building pangenome and pan-model '{label}'...")
    run_pan(
        Path(genomes),
        Path(output),
        label,
        threads=threads,
        add_unique=add_unique,
        soft_core_override=soft_core,
    )


@cli.command()
@click.option(
    "-g",
    "--genomes",
    type=click.Path(exists=True),
    help="Directory of strain genome files (.gbff or .faa format).",
)
@click.option(
    "-f",
    "--fasta",
    type=click.Path(exists=True),
    help="Path to the PPanGGOLiN soft-core FASTA file.",
)
@click.option(
    "-r",
    "--panmodel-rxns",
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
@click.option(
    "--diamond-bits",
    default=150,
    show_default=True,
    help="Minimum bitscore threshold for DIAMOND blastp alignment.",
)
@click.option(
    "--gapseq-bits",
    default=150,
    show_default=True,
    help="Minimum bitscore threshold for gapseq pathway/transport searches.",
)
@click.option(
    "--biomass",
    default="pos",
    show_default=True,
    help="Biomass reaction type for gapseq draft model generation (e.g. 'pos', 'neg').",
)
def strains(
    genomes,
    fasta,
    panmodel_rxns,
    output,
    medium,
    threads,
    diamond_bits,
    gapseq_bits,
    biomass,
):
    """Step 3: Derive strain-specific models using diamond and gapseq."""
    click.echo("Deriving strain models...")

    # Narrow types for static checkers: ensure required args are present
    assert (
        genomes is not None
        and fasta is not None
        and panmodel_rxns is not None
        and medium is not None
    ), "Missing required arguments: genomes, fasta, panmodel_rxns, and medium must be provided"

    derive_strain_models(
        genomes_dir=Path(genomes),
        pangenome_fasta=Path(fasta),
        panmodel_rxn_table=Path(panmodel_rxns),
        output_dir=Path(output),
        medium_file=Path(medium),
        threads=threads,
        diamond_bits=diamond_bits,
        gapseq_bits=gapseq_bits,
        biomass=biomass,
    )


@cli.command()
@click.option(
    "-s",
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
    "-d",
    "--db",
    required=True,
    type=click.Path(exists=True),
    help="Path to the Excel Curation Database containing curated reactions.",
)
def curate(source_model, curated_model, db):
    """Step 4: Apply curation database rules to a specific model."""
    click.echo(f"Curating model {Path(source_model).name}...")
    apply_curations(Path(source_model), Path(db), Path(curated_model))


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


@cli.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-g",
    "--genomes",
    type=click.Path(exists=True),
    help="Directory of new strain genome files (.faa or .gbff format).",
)
@click.option(
    "-f",
    "--fasta",
    type=click.Path(exists=True),
    help="Path to an existing soft-core FASTA file.",
)
@click.option(
    "-r",
    "--panmodel-rxns",
    type=click.Path(exists=True),
    help="Path to an existing pan-model reaction table file.",
)
@click.option(
    "-m",
    "--medium",
    type=click.Path(exists=True),
    help="Path to the media CSV file for gapfilling.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Target output directory for the project workflow.",
)
@click.option(
    "-c",
    "--curation-db",
    type=click.Path(),
    help="Path to the Excel DB for automated curation.",
)
@click.option(
    "-t",
    "--threads",
    default=14,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    "--diamond-bits",
    default=150,
    show_default=True,
    help="Minimum bitscore threshold for DIAMOND blastp alignment.",
)
@click.option(
    "--gapseq-bits",
    default=150,
    show_default=True,
    help="Minimum bitscore threshold for gapseq pathway/transport searches.",
)
@click.option(
    "--biomass",
    default="pos",
    show_default=True,
    help="Biomass reaction type for gapseq draft model generation (e.g. 'pos', 'neg').",
)
@click.option(
    "--no-polish", is_flag=True, help="Skip the ModelPolisher execution step."
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
    "--pause_for_curation",
    is_flag=True,
    help="Pause the pipeline before automatic curation to allow manual evaluation/edits.",
)
@click.option(
    "--resume",
    is_flag=True,
    help="Resume the configured pipeline originating from Step 4 using a finished curation database.",
)
@click.option(
    "-d",
    "--drafts",
    type=click.Path(exists=True, path_type=Path),
    help="Optional path to draft models for zooming/resuming. Defaults to <output>/03-strains/models.",
)
@click.option(
    "--config",
    "config_file",
    type=click.Path(exists=True),
    help="Path to a YAML config file. CLI arguments override config values. Requires pyyaml (pip install 'parmmosahn[config]').",
)
def project(
    genomes,
    fasta,
    panmodel_rxns,
    medium,
    output,
    curation_db,
    threads,
    diamond_bits,
    gapseq_bits,
    biomass,
    no_polish,
    engine,
    pause_for_curation,
    resume,
    drafts,
    config_file,
):
    """
    Rapidly project new strains (.faa or .gbff format) onto an existing pan-model and continue with downstream curation.
    This bypasses Steps 1 & 2 (PPanGGOLiN).
    Can also resume the project using the --resume flag.
    """
    from parmmosahn.curation import generate_curation_template
    from parmmosahn.models import derive_strain_models
    from parmmosahn.pipeline import run_downstream

    if config_file:
        from parmmosahn.utils import load_config

        cfg = load_config(Path(config_file))
        genomes = genomes or cfg.get("genomes")
        fasta = fasta or cfg.get("fasta")
        panmodel_rxns = panmodel_rxns or cfg.get("panmodel_rxns")
        medium = medium or cfg.get("medium")
        threads = threads if threads != 14 else cfg.get("threads", 14)
        diamond_bits = (
            diamond_bits if diamond_bits != 150 else cfg.get("diamond_bits", 150)
        )
        gapseq_bits = gapseq_bits if gapseq_bits != 150 else cfg.get("gapseq_bits", 150)
        biomass = biomass if biomass != "pos" else cfg.get("biomass", "pos")
        engine = engine if engine != "docker" else cfg.get("engine", "docker")

    out_path = Path(output)
    out_path.mkdir(parents=True, exist_ok=True)

    # Handle Resume mode
    if resume:
        if not drafts:
            drafts = out_path / "03-strains" / "models"
            if not drafts.exists():
                raise click.BadParameter(
                    f"Could not find draft models at {drafts}. Please specify --drafts manually."
                )

        if not curation_db:
            curation_db = out_path / "curation_template.xlsx"
        else:
            curation_db = Path(curation_db)

        if not curation_db.exists():
            raise click.FileError(
                str(curation_db),
                hint="Could not find the curation template. If you renamed or moved it, specify its path via --curation-db.",
            )

        click.echo(
            f"Resuming pipeline using defined curation database: {curation_db.name}"
        )

        run_downstream(
            draft_models_dir=drafts,
            output_dir=out_path,
            curation_db=curation_db,
            no_polish=no_polish,
            engine=engine,
        )
        return

    # Normal mode (Not resuming)
    if not all([genomes, fasta, panmodel_rxns, medium]):
        raise click.UsageError(
            "Missing required arguments for normal project strain derivation."
        )

    # 1. Run Strain Derivation (Step 3)
    click.echo(step_banner(3))
    strains_out = out_path / "03-strains"

    # Narrow types for the static type checker: ensure required CLI args are present
    assert (
        genomes is not None
        and fasta is not None
        and panmodel_rxns is not None
        and medium is not None
    ), "Missing required arguments: genomes, fasta, panmodel_rxns, and medium must be provided"

    # Convert to Path objects before passing to the implementation (expects str | Path)
    draft_models_dir = derive_strain_models(
        genomes_dir=Path(genomes),
        pangenome_fasta=Path(fasta),
        panmodel_rxn_table=Path(panmodel_rxns),
        output_dir=strains_out,
        medium_file=Path(medium),
        threads=threads,
        diamond_bits=diamond_bits,
        gapseq_bits=gapseq_bits,
        biomass=biomass,
    )

    if pause_for_curation:
        click.echo("=== Preparing Human-in-the-Loop Curation (running memote) ===")
        generate_curation_template(draft_models_dir, out_path)
        click.echo("=== Pipeline paused for Human-in-the-Loop Curation ===")
        click.echo(
            f"ACTION REQUIRED: Please fill out the curation template at {(out_path / 'curation_template.xlsx').resolve()}"
        )
        click.echo(
            f"When finished, resume using: parmmosahn project --resume -o {out_path.resolve()}"
        )
        return

    # 2. Pipe directly into downstream processing (Steps 4-7)
    run_downstream(
        draft_models_dir=draft_models_dir,
        output_dir=out_path,
        curation_db=Path(curation_db) if curation_db else None,
        no_polish=no_polish,
        engine=engine,
    )


_INIT_CONFIG_TEMPLATE = """\
# PaRMMoSaHN configuration file
# Generated by: parmmosahn init-config
# All values shown are defaults. Remove or comment out any you don't need.
# CLI flags always override values set here.
# You can also use environment variables: PARMMOSAHN_THREADS=32, etc.

# --- Input paths ---
# genomes: /path/to/genomes          # Directory of genome files (.gbff or .faa)
# fasta: /path/to/softcore.faa       # Existing soft-core FASTA (project command)
# panmodel_rxns: /path/to/rxns.tbl   # Existing pan-model reaction table (project command)
# medium: /path/to/medium.csv        # Media CSV file for gapfilling

# --- Output ---
# output: 00-pan_results             # Base output directory

# --- Pangenome (pan command) ---
# label: myproject                   # Label for the pangenome
# add_unique: false                  # Include unique/cloud genes in pan-model
# soft_core: null                    # Override soft-core threshold (fraction 0-1).
#                                    # Default: 2/n_genomes (biologically motivated).

# --- Performance ---
threads: 14
diamond_bits: 150
gapseq_bits: 150
biomass: pos                         # Biomass reaction type: 'pos' or 'neg'

# --- Downstream ---
engine: docker                       # Container engine: docker, podman, or apptainer
"""


@cli.command(name="init-config")
@click.option(
    "-o",
    "--output",
    default="parmmosahn_config.yml",
    show_default=True,
    type=click.Path(),
    help="Path to write the default configuration file.",
)
def init_config(output: str):
    """Generate a default YAML configuration file for reproducible runs."""
    out_path = Path(output)
    if out_path.exists():
        click.confirm(f"{out_path} already exists. Overwrite?", abort=True)
    out_path.write_text(_INIT_CONFIG_TEMPLATE, encoding="utf-8")
    click.echo(f"Configuration template written to {out_path.resolve()}")
    click.echo("Edit the file, then run: parmmosahn project --config " + str(out_path))


@cli.group()
def analyze():
    """Optional analysis extensions: FBA summaries, auxotrophy screening, heatmaps."""
    pass


@analyze.command()
@click.option(
    "-M",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing SBML models (.xml).",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Output TSV file path.",
)
@click.option(
    "-m",
    "--medium",
    type=click.Path(exists=True),
    default=None,
    help="Optional medium CSV (cols: exchange_id, flux) to override model media.",
)
def fba(models, output, medium):
    """Run FBA on each model and produce a comparative growth summary."""
    from parmmosahn.analysis import fba_summary

    click.echo(f"Running FBA summary on models in {models}...")
    fba_summary(
        models_dir=Path(models),
        output_path=Path(output),
        medium_path=Path(medium) if medium else None,
    )
    click.echo(f"FBA summary saved to {output}")


@analyze.command()
@click.option(
    "-M",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing SBML models (.xml).",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Output TSV file path.",
)
@click.option(
    "-m",
    "--medium",
    type=click.Path(exists=True),
    default=None,
    help="Optional medium CSV (cols: exchange_id, flux) to override model media.",
)
@click.option(
    "--threshold",
    default=1e-3,
    show_default=True,
    help="Growth rate below this is flagged as auxotrophic.",
)
def auxotrophy(models, output, medium, threshold):
    """Screen for auxotrophies by single-component medium knockouts."""
    from parmmosahn.analysis import auxotrophy_screen

    click.echo(f"Running auxotrophy screen on models in {models}...")
    auxotrophy_screen(
        models_dir=Path(models),
        output_path=Path(output),
        medium_path=Path(medium) if medium else None,
        threshold=threshold,
    )
    click.echo(f"Auxotrophy screen saved to {output}")


@analyze.command()
@click.option(
    "-M",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing SBML models (.xml).",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Output TSV file path.",
)
@click.option(
    "--plot",
    is_flag=True,
    help="Also generate a clustered heatmap PNG (requires pip install 'parmmosahn[viz]').",
)
def heatmap(models, output, plot):
    """Build a reaction presence/absence matrix across strains."""
    from parmmosahn.analysis import reaction_heatmap

    click.echo(f"Building reaction heatmap from models in {models}...")
    reaction_heatmap(
        models_dir=Path(models),
        output_path=Path(output),
        plot=plot,
    )
    click.echo(f"Reaction matrix saved to {output}")


@cli.command()
@click.option(
    "-M",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing SBML models (.xml) to evaluate.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Directory for memote reports (HTML + JSON).",
)
@click.option(
    "--curation-template",
    is_flag=True,
    help="Also generate the curation Excel template from memote results.",
)
def evaluate(models, output, curation_template):
    """Run memote evaluation on a directory of SBML models."""
    from parmmosahn.curation import (
        extract_memote_problems,
        run_memote_snapshot,
    )

    models_path = Path(models)
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    xml_files = list(models_path.rglob("*.xml"))
    if not xml_files:
        click.echo(f"No XML models found in {models_path}")
        return

    click.echo(f"Evaluating {len(xml_files)} model(s) with memote...")
    for xml in xml_files:
        out_json = output_path / f"{xml.stem}.json"
        out_html = output_path / f"{xml.stem}.html"
        run_memote_snapshot(xml, out_json, out_html)

    if curation_template:
        out_excel = output_path / "curation_template.xlsx"
        extract_memote_problems(output_path, models_path, out_excel)
        click.echo(f"Curation template saved to {out_excel}")

    click.echo(f"Memote reports saved to {output_path}")


@cli.command()
def doctor():
    """Check for required external tools and report their versions."""
    from parmmosahn.utils import check_environment

    REQUIRED = {"ppanggolin", "gapseq", "diamond"}
    ENGINES = {"docker", "podman", "apptainer"}

    click.echo("PaRMMoSaHN Environment Check")
    click.echo("=" * 40)

    results = check_environment()

    click.echo("\nRequired tools:")
    all_required_ok = True
    for name in sorted(REQUIRED):
        info = results[name]
        if info["found"]:
            click.echo(f"  [OK]      {name}: {info['version']}  ({info['path']})")
        else:
            click.echo(f"  [MISSING] {name}: not found in PATH")
            all_required_ok = False

    click.echo("\nContainer engines (at least one needed for ModelPolisher, Step 5):")
    any_engine_ok = False
    for name in sorted(ENGINES):
        info = results[name]
        if info["found"]:
            click.echo(f"  [OK]  {name}: {info['version']}  ({info['path']})")
            any_engine_ok = True
        else:
            click.echo(f"  [--]  {name}: not found")

    click.echo()
    if all_required_ok:
        click.echo("All required tools are available.")
    else:
        click.echo("WARNING: Some required tools are missing. The pipeline will fail.")
    if not any_engine_ok:
        click.echo(
            "WARNING: No container engine found. ModelPolisher (Step 5) will be skipped."
        )


def main():
    cli()


if __name__ == "__main__":
    main()
