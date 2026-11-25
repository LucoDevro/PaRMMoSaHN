import logging
import os
import time
from pathlib import Path

import click

from parmmosahn.curation import apply_curations, init_curation_table
from parmmosahn.models import build_pan_model, derive_strain_models
from parmmosahn.results import gather_and_convert_models
from parmmosahn.pangenome import build_pangenome

logger = logging.getLogger(__name__)


def get_default_threads() -> int:
    """Default to 75 % of available CPU cores (minimum 1)."""
    cpus = os.cpu_count()
    return max(1, int(cpus * 0.75)) if cpus else 1


_DEFAULT_THREADS = get_default_threads()


def _format_elapsed(seconds: float) -> str:
    """Format elapsed seconds into a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s}s"
    h, m = divmod(m, 60)
    return f"{h}h {m}m {s}s"


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
        help_option_names=["-h", "--help"],
        auto_envvar_prefix="PARMMOSAHN",
    )
)
@click.version_option(None, "-v", "--version")
@click.option("--verbose", is_flag=True, help="Show debug-level output.")
@click.option("-q", "--quiet", is_flag=True, help="Only show warnings and errors.")
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
    default=".",
    type=click.Path(),
    show_default=True,
    help="Target base directory for pan-related output (creates 01-pangenome/ and 02-panmodel/ inside).",
)
@click.option(
    "-l",
    "--label",
    required=True,
    help="A descriptive label for the pangenome (e.g., 'clostridiales').",
)
@click.option(
    '-i',
    '--identity',
    default = 0.8,
    help = "Identity threshold for the gene family clustering by PPanGGOLiN",
    show_default = True,
)
@click.option(
    '-c',
    '--coverage',
    default = 0.6,
    help = "Coverage threshold for the gene family clustering by PPanGGOLiN",
    show_default = True,
)
@click.option(
    "-t",
    "--threads",
    default=_DEFAULT_THREADS,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    "--add-small",
    is_flag=True,
    help="Whether to add genes of small gene families (occurrence below threshold) to the pan-model. This can increase coverage but may add noise.",
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
def pan(genomes: Path | None,
        output: Path | None,
        label: str,
        cluster_id: float = 0.80,
        cluster_cov: float = 0.60,
        threads: int = 1,
        add_small: bool = False,
        occurrence: float | None = None,
        ):
    """Step 1 & 2: Construct the partitioned pangenome and build the pan-model."""
    if genomes is None or label is None:
        raise click.UsageError(
            "Missing required arguments: genomes and label must be provided"
        )
        
    output.mkdir(parents=True, exist_ok=True)
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
    )
    
    # --- STEP 2: Pan-model ---
    panmodel_out = output / "02-panmodel"
    build_pan_model(included_fasta, panmodel_out, threads=threads)
    
    logger.info(f"Pan-tasks complete! Pan-model saved to {panmodel_out}")


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
    help="Path to the PPanGGOLiN representatives FASTA file.",
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
    default=_DEFAULT_THREADS,
    show_default=True,
    help="Number of CPU threads to allocate.",
)
@click.option(
    "--no-gapfill",
    is_flag = True,
    help = "Do not gapfill the projected models",
)
@click.option(
    "--strain-bits",
    default=150,
    show_default=True,
    help="Minimum bitscore threshold for DIAMOND blastp alignment.",
)
@click.option(
    "--pan-bits",
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
    "--parallel-strains",
    default=None,
    type=int,
    help="Number of strains to process concurrently. Default: threads//4 (auto). Each worker peaks at ~1.5 GB RAM.",
)
def project(
    genomes,
    fasta,
    panmodel_rxns,
    output,
    medium,
    threads,
    no_gapfill,
    strain_bits,
    pan_bits,
    biomass,
    parallel_strains,
):
    """Step 3: Derive strain-specific models by homology matches."""
    click.echo("Deriving strain models...")

    missing_args = []
    if genomes is None:
        missing_args.append("genomes")
    if fasta is None:
        missing_args.append("fasta")
    if panmodel_rxns is None:
        missing_args.append("panmodel_rxns")
    if medium is None:
        missing_args.append("medium")
    if missing_args:
        raise click.UsageError("Missing required arguments: " + ", ".join(missing_args))

    derive_strain_models(
        genomes_dir=Path(genomes),
        pangenome_fasta=Path(fasta),
        panmodel_rxn_table=Path(panmodel_rxns),
        output_dir=Path(output),
        medium_file=Path(medium),
        threads=threads,
        gapfill = not no_gapfill,
        strain_bits=strain_bits,
        pan_bits=pan_bits,
        biomass=biomass,
        parallel_strains=parallel_strains,
    )
    
    
@cli.command()
@click.option(
    "-o",
    "--output",
    required = True,
    type = click.Path(),
    help = "Path to write the empty curation DB to."
)
def init_curation(output: Path) -> None:
    """
    Write an empty PaRMMoSaHN curation DB.
    """
    output_excel = Path(output)
    logger.info(f"Writing curation template to {output_excel.name}...")
    
    init_curation_table(output_excel)


@cli.command()
@click.option(
    "-s",
    "--source-model",
    required=True,
    type=click.Path(exists=True),
    help="Path to the input draft SBML model (.xml).",
)
@click.option(
    "-o",
    "--output",
    "--curated-model",
    "curated_model",
    required=True,
    type=click.Path(),
    help="Output directory for the curated model (the filename is derived from "
    "the source model). '--curated-model' is a deprecated alias.",
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
    out = Path(curated_model)
    # apply_curations expects an output *directory* and derives the filename
    # from the source model. -o/--output is therefore a directory path.
    out.mkdir(parents=True, exist_ok=True)
    apply_curations(Path(source_model), Path(db), out)


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
    models_path = Path(models)
    xml_count = len(list(models_path.rglob("*.xml")))
    click.echo(f"Input:  {models_path} ({xml_count} models)")
    click.echo(f"Output: {output}")

    t_start = time.time()
    gather_and_convert_models(models, output)
    click.echo(f"Done ({_format_elapsed(time.time() - t_start)} total)")


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
    help="Output TSV file path (additional files written alongside).",
)
@click.option(
    "--plot",
    is_flag=True,
    help="Also generate a dendrogram PNG (requires pip install 'parmmosahn[viz]' scipy).",
)
@click.option(
    "--core-threshold",
    type=float,
    default=0.99,
    show_default=True,
    help="Core if a reaction occurs in strictly more than this fraction of strains (>99% thesis convention).",
)
@click.option(
    "--cloud-threshold",
    type=float,
    default=0.05,
    show_default=True,
    help="Unique (cloud) if a reaction occurs in strictly fewer than this fraction of strains (<5%).",
)
def panreactome(models, output, plot, core_threshold, cloud_threshold):
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


@analyze.command(name="prune-report")
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
def prune_report(models, output):
    """Report dead-end metabolites + fully-blocked reactions per model.

    This is a DRY-RUN diagnostic — no SBML files are modified. It identifies
    metabolites with zero producers or zero consumers (the structural source
    of the cohort-wide ~42% blocked-fraction baseline) and reactions whose
    participants are all dead-end (guaranteed-prunable without affecting any
    connected pathway).

    Scaffolding for the full `parmmosahn prune` step on the v0.3 roadmap
    (item B2); use the output to size the prune-step impact before
    committing to a destructive pass.
    """
    from parmmosahn.analysis import prune_report as prune_report_fn

    click.echo(f"Building prune dry-run report from models in {models}...")
    prune_report_fn(models_dir=Path(models), output_path=Path(output))
    click.echo(f"Prune dry-run report saved to {output}")


@analyze.command()
@click.option(
    "-M",
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
    "--threshold",
    default=1e-3,
    show_default=True,
    help="Growth rate below this is predicted as 'no growth'.",
)
def validate(models, phenotypes, output, medium, threshold):
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
    "--essential-fraction",
    default=1e-2,
    show_default=True,
    help="Knockout growth below this fraction of WT marks a gene as essential.",
)
@click.option(
    "--min-wt-growth",
    default=1e-3,
    show_default=True,
    help="Absolute wild-type growth floor below which a strain is skipped.",
)
def essentiality(models, output, medium, essential_fraction, min_wt_growth):
    """Predict essential genes via single-gene deletion FBA."""
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
@click.option(
    "--workers",
    default=1,
    show_default=True,
    type=int,
    help="Number of strains to evaluate in parallel. Memote uses ~6 internal "
    "threads per model, so 2-3 outer workers typically saturate a 12-16 core host.",
)
@click.option(
    "--skip-existing/--no-skip-existing",
    default=True,
    show_default=True,
    help="Skip strains whose output JSON already exists and is non-empty. "
    "Allows safe restart after a crash without redoing completed work.",
)
def evaluate(models, output, curation_template, workers, skip_existing):
    """Run memote evaluation on a directory of SBML models."""
    from parmmosahn.curation import (
        extract_memote_problems,
        run_memote_snapshot,
    )

    models_path = Path(models)
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    all_xml = sorted(models_path.rglob("*.xml"))
    if not all_xml:
        click.echo(f"No XML models found in {models_path}")
        return

    # Outputs are keyed by stem ({stem}.json/.html), so two .xml files sharing
    # a stem would both target the same report and be evaluated as separate
    # tasks. This happens with gapseq trees, which keep a flat
    # {strain}.xml *and* a per-strain working subdir containing more SBML:
    # rglob then returns the stem 2-3x. --skip-existing cannot dedupe these
    # (neither output exists yet when both are queued), so with --workers>1
    # the same strain runs concurrently, twice. Collapse by stem here, first
    # path wins, so each strain is evaluated exactly once.
    xml_files = []
    seen_stems: set[str] = set()
    dupe_stems: set[str] = set()
    for xml in all_xml:
        if xml.stem in seen_stems:
            dupe_stems.add(xml.stem)
            continue
        seen_stems.add(xml.stem)
        xml_files.append(xml)

    n_total = len(xml_files)
    if dupe_stems:
        click.echo(
            f"Note: {len(all_xml)} .xml files collapse to {n_total} unique "
            f"strains ({len(dupe_stems)} stem(s) had extra copies in "
            f"sub-directories; first path used)."
        )

    todo = []
    for xml in xml_files:
        out_json = output_path / f"{xml.stem}.json"
        if skip_existing and out_json.exists() and out_json.stat().st_size > 0:
            continue
        todo.append(xml)

    n = len(todo)
    skipped = n_total - n
    click.echo(f"Input:  {models_path} ({n_total} models)")
    click.echo(f"Output: {output_path}")
    if skipped:
        click.echo(f"Skipping {skipped} strain(s) with existing JSON output.")

    if n == 0:
        click.echo("Nothing to do.")
        if curation_template:
            out_excel = output_path / "curation_template.xlsx"
            extract_memote_problems(output_path, models_path, out_excel)
            click.echo(f"Curation template saved to {out_excel}")
        return

    # Sanity-check RAM before launching parallel memote workers. memote peaks
    # around 3-4 GB per worker on typical strain models; on memory-constrained
    # hosts (e.g. default WSL2 caps) requesting too many workers triggers
    # BrokenProcessPool errors. Suggest a safer count instead of letting the
    # pool fail mid-run. See suggestions_code.md §5.
    from parmmosahn.utils import get_system_resources

    resources = get_system_resources()
    avail_gb = resources.get("available_ram_gb")
    if avail_gb is not None and workers > 1:
        safe_workers = max(1, int(avail_gb // 4))
        if workers > safe_workers:
            click.echo(
                f"WARNING: requested workers={workers} but only "
                f"{avail_gb:.1f} GB RAM available "
                f"(~4 GB/worker recommended). Consider `--workers {safe_workers}` "
                f"if memote crashes with BrokenProcessPool."
            )

    click.echo(f"Evaluating {n} model(s) with memote (workers={workers})...")
    t_start = time.time()

    failed = 0
    if workers <= 1:
        for i, xml in enumerate(todo, 1):
            t_model = time.time()
            click.echo(f"  [{i}/{n}] {xml.stem} ...", nl=False)
            out_json = output_path / f"{xml.stem}.json"
            out_html = output_path / f"{xml.stem}.html"
            ok = run_memote_snapshot(xml, out_json, out_html)
            if not ok:
                failed += 1
            click.echo(
                f" {'done' if ok else 'FAIL'} "
                f"({_format_elapsed(time.time() - t_model)})"
            )
    else:
        # Threads orchestrate; each task delegates to run_memote_snapshot which
        # spawns its own subprocess (subprocess_mode=True by default), so memory
        # is reclaimed per strain even though the outer pool is thread-based.
        from concurrent.futures import ThreadPoolExecutor, as_completed

        futures = {}
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for xml in todo:
                out_json = output_path / f"{xml.stem}.json"
                out_html = output_path / f"{xml.stem}.html"
                fut = executor.submit(run_memote_snapshot, xml, out_json, out_html)
                futures[fut] = xml
            for done, fut in enumerate(as_completed(futures), 1):
                xml = futures[fut]
                t_done = time.time() - t_start
                try:
                    ok = fut.result()
                    status = "ok" if ok else "FAIL"
                except Exception as e:
                    ok = False
                    status = f"ERROR: {e}"
                if not ok:
                    failed += 1
                click.echo(
                    f"  [{done}/{n}] {xml.stem} {status} "
                    f"(elapsed {_format_elapsed(t_done)})"
                )

    if curation_template:
        out_excel = output_path / "curation_template.xlsx"
        extract_memote_problems(output_path, models_path, out_excel)
        click.echo(f"Curation template saved to {out_excel}")

    click.echo(
        f"Memote reports saved to {output_path} ({_format_elapsed(time.time() - t_start)} total)"
    )
    if failed:
        # Signal failure so CI gates / `&&` chains can detect that strains failed.
        click.echo(f"{failed}/{n} model(s) failed memote evaluation.")
        raise SystemExit(1)


@cli.command()
@click.option(
    "-M",
    "--models",
    required=True,
    type=click.Path(exists=True),
    help="Directory containing SBML models (.xml) to enrich.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(),
    help="Directory for the annotation-enriched models.",
)
@click.option(
    "--rename-to-bigg",
    is_flag=True,
    default=False,
    help=(
        "After annotation enrichment, rewrite SEED-namespace primary IDs to "
        "BiGG IDs where a mapping exists (metabolites, exchanges, internal "
        "reactions). Opt-in; default off. Unmapped IDs keep their SEED form."
    ),
)
def polish(models, output, rename_to_bigg):
    """Enrich model cross-ref annotations."""
    from parmmosahn.polish import enrich_models_directory

    models_path = Path(models)
    output_path = Path(output)

    xml_count = len(list(models_path.rglob("*.xml")))
    click.echo(f"Input:  {models_path} ({xml_count} models)")
    click.echo(f"Output: {output_path}")
    if rename_to_bigg:
        click.echo("Rename: SEED -> BiGG primary IDs (--rename-to-bigg)")

    t_start = time.time()
    enrich_models_directory(models_path, output_path, rename=rename_to_bigg)
    click.echo(
        f"Enriched models saved to {output_path} ({_format_elapsed(time.time() - t_start)} total)"
    )


_DOCTOR_REQUIRED_TOOLS = ("ppanggolin", "gapseq", "diamond")
_DOCTOR_ENGINES = ("docker", "podman", "apptainer")


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
    engines = {name: env[name] for name in _DOCTOR_ENGINES}

    avail_gb = resources.get("available_ram_gb")
    memote_workers_recommended = (
        max(1, int(avail_gb // 4)) if avail_gb is not None else None
    )

    all_required_ok = all(info["found"] for info in required.values())
    any_engine_ok = any(info["found"] for info in engines.values())
    all_data_ok = all(info["found"] for info in bundled.values())

    warnings: list[str] = []
    if not all_required_ok:
        warnings.append(
            "Some required external tools are missing. The pipeline will fail."
        )
    if not any_engine_ok:
        warnings.append(
            "No container engine found. ModelPolisher (Step 5b) will be skipped."
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
        "container_engines": engines,
        "bundled_data": bundled,
        "host_resources": resources,
        "memote_workers_recommended": memote_workers_recommended,
        "checks": {
            "all_required_tools_ok": all_required_ok,
            "any_container_engine_ok": any_engine_ok,
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

    click.echo("\nContainer engines (at least one needed for ModelPolisher, Step 5b):")
    for name in sorted(report["container_engines"]):
        info = report["container_engines"][name]
        if info["found"]:
            click.echo(f"  [OK]  {name}: {info['version']}  ({info['path']})")
        else:
            click.echo(f"  [--]  {name}: not found")

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


@cli.command()
@click.option(
    "--json",
    "json_out",
    is_flag=True,
    default=False,
    help=(
        "Emit the environment report as a single JSON object on stdout "
        "instead of the human-readable text. Useful for CI gates and "
        "automation (pipe to `jq` to filter on specific tools / versions)."
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
