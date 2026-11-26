import logging
import time
import shutil
from pathlib import Path

import click

from parmmosahn.pangenome import build_pangenome
from parmmosahn.models import build_pan_model, derive_strain_models
from parmmosahn.curation import evaluate_models, extract_memote_problems, apply_curations, init_curation_table
from parmmosahn.polish import polish_models
from parmmosahn.convert import convert_models
from parmmosahn.utils import format_elapsed, safe_init_output_dir, get_default_threads, SpecialHelpOrder

logger = logging.getLogger(__name__)


def _setup_logging(verbose: bool, quiet: bool) -> None:
    """Configure logging level based on CLI flags."""
    if verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s [%(name)s] %(message)s",
        datefmt = "%H:%M:%S",
        force=True,
    )

    # Suppress noisy cobra solver warnings that clutter output
    logging.getLogger("cobra.core.model").setLevel(logging.WARNING)


@click.group(
    context_settings=dict(
        max_content_width = shutil.get_terminal_size().columns - 10,
        help_option_names=["-h", "--help"],
        auto_envvar_prefix="PARMMOSAHN",
    ),
    cls = SpecialHelpOrder,
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Only show warnings and errors."
)
@click.option(
    "--verbose",
    is_flag=True,
    help="Show debug-level output."
)
@click.version_option(None, "-v", "--version")
@click.pass_context
def cli(ctx, verbose, quiet):
    """
    PaRMMoSaHN: Pangenome Reference-based Metabolic Modelling by Saving Homolog Networks.

    All options can also be set via environment variables prefixed with PARMMOSAHN_.
    For example: PARMMOSAHN_THREADS=32, PARMMOSAHN_BIOMASS=neg.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet
    _setup_logging(verbose, quiet)


@cli.command(help_priority = 1)
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
    default="output",
    type=click.Path(),
    show_default=True,
    help="Target base directory for pan-related output (creates 01-pangenome/ and 02-panmodel/ inside).",
)
@click.option(
    "-l",
    "--label",
    default='reference',
    show_default=True,
    help="A descriptive label for the pangenome.",
)
@click.option(
    '-i',
    '--identity',
    'cluster_id',
    default = 0.8,
    help = "Identity threshold for the gene family clustering by PPanGGOLiN",
    show_default = True,
)
@click.option(
    '-c',
    '--coverage',
    'cluster_cov',
    default = 0.6,
    help = "Coverage threshold for the gene family clustering by PPanGGOLiN",
    show_default = True,
)
@click.option(
    "--add-small",
    is_flag=True,
    help="Add genes of small gene families (with occurrence below threshold) to the panmodel. This can increase reaction coverage but may add noise.",
)
@click.option(
    "--occurrence",
    default=None,
    type=float,
    show_default=True,
    help=(
        "Override the default gene family occurrence threshold (fraction 0–1). "
        "Default is 2/n_genomes."
    ),
)
@click.option(
    "-t",
    "--threads",
    default=get_default_threads(),
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def pan(genomes: Path | None,
        output: Path | None,
        label: str,
        cluster_id: float = 0.80,
        cluster_cov: float = 0.60,
        threads: int = 1,
        add_small: bool = False,
        occurrence: float | None = None,
        force: bool = False,
        ):
    """Construct the pangenome and build the panmodel."""
    if genomes is None or label is None:
        raise click.UsageError(
            "Missing required arguments: genomes and label must be provided"
        )
    
    safe_init_output_dir(output, force = force)
    click.echo(f"Building pangenome and pan-model '{label}'...")
    
    # --- STEP 1: Pangenome ---
    pangenome_out = output / "01-pangenome"
    included_fasta = build_pangenome(
        genomes,
        pangenome_out,
        label,
        cluster_id = cluster_id,
        cluster_cov = cluster_cov,
        threads=threads,
        add_small=add_small,
        occurrence=occurrence,
        force = force,
    )
    
    # --- STEP 2: Pan-model ---
    panmodel_out = output / "02-panmodel"
    build_pan_model(included_fasta, panmodel_out, threads=threads, force = force)
    
    logger.info(f"Pan-tasks complete! Pan-model saved to {panmodel_out}")


@cli.command(help_priority = 2)
@click.option(
    "-sg",
    "--strains",
    type=click.Path(exists=True),
    help="Directory of strain genome files (.gbff or .faa format).",
)
@click.option(
    "-rf",
    "--reference",
    type=click.Path(exists=True),
    help="Path to the PPanGGOLiN representatives FASTA file.",
)
@click.option(
    "-pr",
    "--panmodel",
    'panmodel_rxns',
    type=click.Path(exists=True),
    help="Path to the panmodel reaction table.",
)
@click.option(
    "-o",
    "--output",
    default="output/03-strain_models",
    type=click.Path(),
    show_default=True,
    help="Target directory for the generated strain-specific draft models.",
)
@click.option(
    "-m",
    "--medium",
    type=click.Path(exists=True),
    help="Medium CSV file used for gapfilling.",
)
@click.option(
    "--strain-bits",
    default=150,
    show_default=True,
    help="Minimum bitscore threshold for DIAMOND homolog detection in the projected strains.",
)
@click.option(
    "--pan-bits",
    default=150,
    show_default=True,
    help="Minimum bitscore threshold for gapseq panmodel pathway/transport identification.",
)
@click.option(
    "--biomass",
    default="pos",
    type = click.Choice(['pos', 'neg']),
    show_default=True,
    show_choices=True,
    help="Biomass type for gapseq draft model generation ('pos' or 'neg').",
)
@click.option(
    "--no-gapfill",
    is_flag = True,
    help = "Do not gapfill the projected models",
)
@click.option(
    "-t",
    "--threads",
    default=get_default_threads(),
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    '-w',
    "--workers",
    default=None,
    type=int,
    help="Number of concurrent workers. Default: threads//4. Account for ~1.5 GB RAM per worker.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def project(
    strains,
    reference,
    panmodel_rxns,
    output,
    medium,
    threads,
    no_gapfill,
    strain_bits,
    pan_bits,
    biomass,
    workers,
    force
):
    """Project strain-specific models from the panmodel by homology matches."""
    click.echo("Deriving strain models...")

    missing_args = []
    if strains is None:
        missing_args.append("genomes")
    if reference is None:
        missing_args.append("fasta")
    if panmodel_rxns is None:
        missing_args.append("panmodel_rxns")
    if medium is None:
        missing_args.append("medium")
    if missing_args:
        raise click.UsageError("Missing required arguments: " + ", ".join(missing_args))

    derive_strain_models(
        genomes_dir=Path(strains),
        pangenome_fasta=Path(reference),
        panmodel_rxn_table=Path(panmodel_rxns),
        output_dir=Path(output),
        medium_file=Path(medium),
        threads=threads,
        gapfill = not no_gapfill,
        strain_bits=strain_bits,
        pan_bits=pan_bits,
        biomass=biomass,
        parallel_strains=workers,
        force = force,
    )
    
    
@cli.command(help_priority = 3)
@click.option(
    "-i",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Input directory containing SBML models (.xml) to evaluate.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Output directory for memote reports (HTML + JSON).",
)
@click.option(
    '-e',
    "--excel",
    default = None,
    show_default=True,
    help="Path to optional curation Excel prefilled with memote results.",
)
@click.option(
    '-w',
    "--workers",
    default=1,
    show_default=True,
    type=int,
    help="Number of concurrent workers. Memote uses ~6 threads per worker internally.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def evaluate(models, output, excel, workers, force):
    """Run memote evaluation for a directory of models."""
    models_path = Path(models)
    output_path = Path(output)
    
    t_start = time.time()
    failed = evaluate_models(models, output, workers, force)
    click.echo(
        f"Memote reports saved to {output_path} ({format_elapsed(time.time() - t_start)} total)"
    )
    
    if failed:
        click.echo(f"{failed} model(s) failed memote evaluation.")

    if excel:
        extract_memote_problems(output_path, models_path, excel)
        click.echo(f"Curation template saved to {excel}")
    

@cli.command(help_priority = 4)
@click.option(
    "-o",
    "--output",
    required = True,
    type = click.Path(),
    help = "Output path for the empty curation DB."
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def init_curation(output: Path, force: bool) -> None:
    """
    Write an empty curation spreadsheet.
    """
    output_excel = Path(output)
    logger.info(f"Writing curation template to {output_excel.name}...")
    
    init_curation_table(output_excel, force = force)


@cli.command(help_priority = 5)
@click.option(
    "-i",
    "--input",
    'source_models',
    required=True,
    type=click.Path(exists=True),
    help="Path to the input draft SBML model (.xml).",
)
@click.option(
    "-o",
    "--output",
    "curated_models",
    required=True,
    type=click.Path(),
    help="Output directory for the curated model (the filename is derived from "
    "the source model).",
)
@click.option(
    "-e",
    "--excel",
    'db',
    required=True,
    type=click.Path(exists=True),
    help="Excel curation database with manual curations.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def curate(source_models, curated_models, db, force):
    """Curate a directory of models using an existing curation spreadsheet."""
    out = Path(curated_models)
    apply_curations(Path(source_models), Path(db), out, force = force)


@cli.command(help_priority = 6)
@click.option(
    "-i",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing SBML models (.xml) to polish.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Directory for the polished models.",
)
@click.option(
    "--bigg",
    'rename_to_bigg',
    is_flag=True,
    show_default=False,
    help=(
        "After polishing, translate SEED IDs to "
        "BiGG IDs where possible. Unmapped IDs keep their SEED form."
    ),
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def polish(models, output, rename_to_bigg, force):
    """Polish a folder of models by enriching cross-refs."""
    models_path = Path(models)
    output_path = Path(output)

    xml_count = len(list(models_path.rglob("*.xml")))
    click.echo(f"Input:  {models_path} ({xml_count} models)")
    click.echo(f"Output: {output_path}")
    if rename_to_bigg:
        click.echo("Rename: SEED -> BiGG primary IDs (--rename-to-bigg)")

    t_start = time.time()
    polish_models(models_path, output_path, rename=rename_to_bigg, force = force)
    click.echo(
        f"Enriched models saved to {output_path} ({format_elapsed(time.time() - t_start)} total)"
    )


@cli.command(help_priority = 7)
@click.option(
    "-i",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Folder containing SBML models.",
)
@click.option(
    "-o",
    "--output",
    default="output/models",
    type=click.Path(),
    show_default=True,
    help="Target folder for the new output formats (JSON/MATLAB).",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def convert(models, output, force):
    """Convert models to JSON and MATLAB formats."""
    models_path = Path(models)
    xml_count = len(list(models_path.rglob("*.xml")))
    click.echo(f"Input:  {models_path} ({xml_count} models)")
    click.echo(f"Output: {output}")

    t_start = time.time()
    convert_models(models, output, force)
    click.echo(f"Done ({format_elapsed(time.time() - t_start)} total)")


@cli.group(help_priority = 8)
def analyze():
    """Optional analysis extensions: FBA summaries, auxotrophy screening, heatmaps."""
    pass


@analyze.command()
@click.option(
    "-i",
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
    help="Optional medium CSV to override model media.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def fba(models, output, medium, force):
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
    "-i",
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
    "--auxo",
    'threshold',
    default=1e-3,
    show_default=True,
    help="A growth rate below this threshold is flagged as auxotrophic.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def auxotrophy(models, output, medium, threshold, force):
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
    "-i",
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
    '-p',
    "--plot",
    default = None,
    help="Optional clustered heatmap output path (image format deduced from extension).",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def heatmap(models, output, plot, force):
    """Build a reaction presence/absence matrix across strains."""
    from parmmosahn.analysis import reaction_heatmap

    click.echo(f"Building reaction heatmap from models in {models}...")
    reaction_heatmap(
        models_dir=Path(models),
        output_path=Path(output),
        plot=plot,
    )
    click.echo(f"Reaction matrix saved to {output}")


#TODO: to be revisited
@analyze.command()
@click.option(
    "-i",
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
    help="Output TSV file path (additional files written alongside).",
)
@click.option(
    '-p',
    "--plot",
    is_flag=True,
    help="Also generate a dendrogram PNG (requires pip install 'parmmosahn[viz]' scipy).",
)
@click.option(
    "--core",
    'core_threshold',
    type=float,
    default=0.99,
    show_default=True,
    help="Core if a reaction occurs in strictly more than this fraction of strains (>99% thesis convention).",
)
@click.option(
    "--cloud",
    'cloud_threshold',
    type=float,
    default=0.05,
    show_default=True,
    help="Unique (cloud) if a reaction occurs in strictly fewer than this fraction of strains (<5%).",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def panreactome(models, output, plot, core_threshold, cloud_threshold, force):
    """Characterize the pan-reactome: core/accessory/unique reactions, Jaccard similarity, accumulation curve."""
    from parmmosahn.analysis import panreactome as panreactome_fn

    click.echo(f"Analyzing pan-reactome from models in {models}...")
    panreactome_fn(
        models_dir=Path(models),
        output_path=Path(output),
        plot=plot,
        core_threshold=core_threshold,
        cloud_threshold=cloud_threshold,
    )
    click.echo(f"Pan-reactome analysis saved alongside {output}")


@analyze.command(name="dead-ends")
@click.option(
    "-i",
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
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def dead_ends(models, output, force):
    """Report dead-end metabolites and fully-blocked reactions in each model."""
    from parmmosahn.analysis import prune_report as prune_report_fn

    click.echo(f"Building dead-ends report from models in {models}...")
    prune_report_fn(models_dir=Path(models), output_path=Path(output))
    click.echo(f"Dead-ends report saved to {output}")


@analyze.command()
@click.option(
    "-i",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing SBML models (.xml).",
)
@click.option(
    "-p",
    "--phenotypes",
    required=True,
    type=click.Path(exists=True),
    help="CSV/TSV with experimental growth phenotypes (matrix or long format).",
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
    help="Optional base medium CSV (cols: exchange_id, flux).",
)
@click.option(
    '-mg',
    "--min-growth",
    'threshold',
    default=1e-3,
    show_default=True,
    help="Growth rate below this is considered as 'no growth'.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def validate(models, phenotypes, output, medium, threshold, force):
    """Validate model predictions against experimental growth phenotypes."""
    from parmmosahn.analysis import validate_phenotypes

    click.echo(f"Validating models in {models} against {phenotypes}...")
    validate_phenotypes(
        models_dir=Path(models),
        phenotype_path=Path(phenotypes),
        output_path=Path(output),
        medium_path=Path(medium) if medium else None,
        threshold=threshold,
    )
    click.echo(f"Validation results saved to {output}")


@analyze.command()
@click.option(
    "-i",
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
    '-mf',
    "--min-frac",
    'essential_fraction',
    default=1e-2,
    show_default=True,
    help="Minimum fraction of wildtype growth rate for a knockout to be considered as non-essential.",
)
@click.option(
    '-mwr',
    "--min-wt-rate",
    'min_wt_growth',
    default=1e-3,
    show_default=True,
    help="Minimum wildtype growth rate for a strain to be included in the analysis.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    help = "Overwrite existing output"
)
def essentiality(models, output, medium, essential_fraction, min_wt_growth, force):
    """Predict essential genes via single knockout FBA."""
    from parmmosahn.analysis import gene_essentiality

    click.echo(f"Running gene essentiality analysis on models in {models}...")
    gene_essentiality(
        models_dir=Path(models),
        output_path=Path(output),
        medium_path=Path(medium) if medium else None,
        essential_fraction=essential_fraction,
        min_wt_growth=min_wt_growth,
    )
    click.echo(f"Gene essentiality predictions saved to {output}")


_DOCTOR_REQUIRED_TOOLS = ("ppanggolin", "gapseq", "diamond")


def _collect_doctor_report() -> dict:
    """Single source of truth for the data both --text and --json render.

    Splitting the data collection from the rendering keeps the two output
    modes from drifting and makes the JSON shape testable independently of
    the human-readable layout.
    """
    from parmmosahn.utils import (
        check_bundled_data,
        check_environment,
        get_system_resources,
    )

    env = check_environment()
    bundled = check_bundled_data()
    resources = get_system_resources()

    required = {name: env[name] for name in _DOCTOR_REQUIRED_TOOLS}

    avail_gb = resources.get("available_ram_gb")
    memote_workers_recommended = (
        max(1, int(avail_gb // 4)) if avail_gb is not None else None
    )

    all_required_ok = all(info["found"] for info in required.values())
    all_data_ok = all(info["found"] for info in bundled.values())

    warnings: list[str] = []
    if not all_required_ok:
        warnings.append(
            "Some required external tools are missing. The pipeline will fail."
        )
    if avail_gb is not None and avail_gb < 8:
        warnings.append(
            f"Only {avail_gb:.1f} GB RAM available — memote evaluation may "
            f"OOM. Run with `--workers 1` and consider setting "
            f"MALLOC_ARENA_MAX=2 in the environment."
        )
    if not all_data_ok:
        warnings.append(
            "One or more bundled TSV lookup tables are missing from the "
            "installed package. Annotation enrichment (Step 5) will fail. "
            "Re-install with `pip install --force-reinstall -e .` "
            "or rebuild the wheel."
        )

    from importlib.metadata import version as _pkg_version

    return {
        "parmmosahn_version": _pkg_version("parmmosahn"),
        "required_tools": required,
        "bundled_data": bundled,
        "host_resources": resources,
        "memote_workers_recommended": memote_workers_recommended,
        "checks": {
            "all_required_tools_ok": all_required_ok,
            "all_bundled_data_ok": all_data_ok,
        },
        "warnings": warnings,
    }


def _render_doctor_text(report: dict) -> None:
    """Print the human-readable doctor report."""
    click.echo("PaRMMoSaHN Environment Check")
    click.echo("=" * 40)

    click.echo("\nRequired tools:")
    for name in sorted(report["required_tools"]):
        info = report["required_tools"][name]
        if info["found"]:
            click.echo(f"  [OK]      {name}: {info['version']}  ({info['path']})")
        else:
            click.echo(f"  [MISSING] {name}: not found in PATH")

    click.echo("\nBundled data files (annotation enrichment):")
    for name in sorted(report["bundled_data"]):
        info = report["bundled_data"][name]
        if info["found"]:
            kb = info["size_bytes"] / 1024 if info["size_bytes"] else 0
            click.echo(f"  [OK]      {name}: {info['rows']} rows, {kb:.1f} KB")
        else:
            click.echo(f"  [MISSING] {name}: not in the installed package")

    resources = report["host_resources"]
    click.echo("\nHost resources:")
    if "cpu_count" in resources:
        click.echo(f"  CPUs:           {resources['cpu_count']}")
    if "total_ram_gb" in resources:
        click.echo(f"  Total RAM:      {resources['total_ram_gb']:.1f} GB")
    if "available_ram_gb" in resources:
        click.echo(f"  Available RAM:  {resources['available_ram_gb']:.1f} GB")
    if "platform" in resources:
        click.echo(f"  Platform:       {resources['platform']}")
    safe_workers = report["memote_workers_recommended"]
    if safe_workers is not None:
        click.echo(
            f"  Memote workers: recommended max {safe_workers} "
            f"(~4 GB per worker; lower on systems where other jobs share RAM)"
        )

    click.echo()
    if report["checks"]["all_required_tools_ok"]:
        click.echo("All required tools are available.")
    for w in report["warnings"]:
        click.echo(f"WARNING: {w}")


@cli.command(help_priority = 9)
@click.option(
    "--json",
    "json_out",
    is_flag=True,
    default=False,
    help=(
        "Emit the environment report as a single JSON object on stdout "
        "instead of the human-readable text."
    ),
)
def doctor(json_out: bool):
    """Check for required external tools and report their versions."""
    import json as _json

    report = _collect_doctor_report()
    if json_out:
        click.echo(_json.dumps(report, indent=2, sort_keys=True, default=str))
    else:
        _render_doctor_text(report)


def main():
    cli()


if __name__ == "__main__":
    main()
