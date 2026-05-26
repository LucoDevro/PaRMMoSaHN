import logging
import os
import time
from pathlib import Path

import click

from parmmosahn.constants import STEPS, step_banner
from parmmosahn.curation import apply_curations
from parmmosahn.models import build_pan_model, derive_strain_models
from parmmosahn.pangenome import build_pangenome
from parmmosahn.results import gather_and_convert_models

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


def _downstream_step_entry(result: object, elapsed: float) -> dict:
    """Build the ``04-07 downstream`` summary entry from run_downstream's result.

    Records an honest status (``completed_with_failures`` when any strain failed
    curation or fewer final models were produced than expected) plus the polish
    outcome. Tolerant of a non-dict ``result`` (e.g. a mocked run_downstream),
    in which case it degrades to the previous plain ``completed`` entry.
    """
    ds = result if isinstance(result, dict) else {}
    has_failures = bool(ds.get("curation_failed")) or (
        ds.get("n_final", 0) < ds.get("n_expected", 0)
    )
    entry: dict = {
        "status": "completed_with_failures" if has_failures else "completed",
        "elapsed": _format_elapsed(elapsed),
    }
    if ds.get("polish"):
        entry["polish"] = ds["polish"]
    if ds.get("curation_failed"):
        entry["curation_failed_count"] = len(ds["curation_failed"])
    if "n_final" in ds and "n_expected" in ds:
        entry["models_final"] = ds["n_final"]
        entry["models_expected"] = ds["n_expected"]
    return entry


def _stdin_is_tty() -> bool:
    """Indirection so tests can force-toggle the interactive-pause behaviour."""
    import sys

    return sys.stdin.isatty()


def _echo_and_log(message: str) -> None:
    """Echo a message to stdout AND tee it into any active run.log file.

    Plain ``click.echo`` writes only to stdout, which means user-facing
    progress lines emitted from ``run`` / ``project`` were missing from
    ``run.log`` even though logger.info() messages from pipeline modules
    were captured. This helper keeps the click.echo terminal output and
    additionally emits the same message via every FileHandler attached
    to the root logger — but NOT via the basicConfig StreamHandler, so
    nothing double-prints on the terminal.

    Used for status lines emitted from inside ``run`` / ``project``.
    Other commands have no run.log, so the FileHandler loop is a no-op
    and behaviour matches plain ``click.echo``.
    """
    click.echo(message)
    record = logging.LogRecord(
        name="parmmosahn.cli",
        level=logging.INFO,
        pathname="",
        lineno=0,
        msg=message,
        args=None,
        exc_info=None,
    )
    for handler in logging.getLogger().handlers:
        if isinstance(handler, logging.FileHandler):
            handler.handle(record)


def _setup_logging(verbose: bool, quiet: bool) -> None:
    """Configure logging level based on CLI flags."""
    if verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    else:
        level = logging.INFO

    # ``force=True`` re-attaches the StreamHandler to the *current* sys.stderr
    # on every invocation. Without it, a handler bound to a stale stderr from
    # an earlier ``CliRunner.invoke`` (or a long-running interpreter session)
    # would silently swallow logger output. This is what made the v0.2.1
    # ``--pause-for-curation`` notices appear to vanish under ``conda run``.
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
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
    "-t",
    "--threads",
    default=_DEFAULT_THREADS,
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
        # --output defaults to ".", which is always truthy, so a plain `or`
        # would always win and silently ignore a config-supplied output. Apply
        # the config value only when --output was not given on the command line.
        src = click.get_current_context().get_parameter_source("output")
        if src is None or src.name == "DEFAULT":
            output = cfg.get("output", output)
        label = label or cfg.get("label")
        threads = (
            threads
            if threads != _DEFAULT_THREADS
            else cfg.get("threads", _DEFAULT_THREADS)
        )
        add_unique = add_unique or cfg.get("add_unique", False)
        if soft_core is None:
            soft_core = cfg.get("soft_core")
    # Ensure the required values are available (may be provided via config)
    if genomes is None or label is None:
        raise click.UsageError(
            "Missing required arguments: genomes and label must be provided either via CLI or config"
        )

    click.echo(f"Building pangenome and pan-model '{label}'...")
    run_pan(
        Path(genomes),
        Path(output),
        label,
        threads=threads,
        add_unique=add_unique,
        soft_core_override=soft_core,
    )


@cli.command(name="run")
@click.option(
    "-g",
    "--genomes",
    type=click.Path(exists=True),
    help="Directory of input genome files (.gbff format). Required (CLI or config).",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    help="Target base directory for all pipeline output. Required (CLI or config).",
)
@click.option(
    "-l",
    "--label",
    help="A descriptive label for the pangenome (e.g., 'clostridiales'). Required.",
)
@click.option(
    "-m",
    "--medium",
    type=click.Path(exists=True),
    help="Path to the media CSV file used for gapfilling. Required (CLI or config).",
)
@click.option(
    "-t",
    "--threads",
    default=_DEFAULT_THREADS,
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
    "--add-unique",
    is_flag=True,
    help="Include unique genes (not in soft-core) in the pan-model.",
)
@click.option(
    "--soft-core",
    default=None,
    type=float,
    help="Override the soft-core frequency threshold (fraction 0-1).",
)
@click.option(
    "--polish",
    is_flag=True,
    default=False,
    help="Enable ModelPolisher (requires Docker/Podman/Apptainer). Off by default.",
)
@click.option(
    "-e",
    "--engine",
    type=click.Choice(["docker", "podman", "apptainer"]),
    default="docker",
    show_default=True,
    help="Container engine for ModelPolisher (only relevant with --polish).",
)
@click.option(
    "-c",
    "--curation-db",
    type=click.Path(exists=True),
    help="Path to an existing Excel curation database for automated curation.",
)
@click.option(
    "--pause-for-curation",
    is_flag=True,
    help="Pause after memote evaluation to allow manual curation edits before continuing.",
)
@click.option(
    "--parallel-strains",
    default=None,
    type=int,
    help=(
        "Number of strains to process concurrently in Step 3. "
        "Default: threads//4 (auto). Each worker gets threads//N CPU threads "
        "and peaks at ~1.5 GB RAM."
    ),
)
@click.option(
    "--config",
    "config_file",
    type=click.Path(exists=True),
    help=(
        "Path to a YAML config file (see `parmmosahn init-config`). CLI flags "
        "override config values; required args (genomes/output/label/medium) "
        "may be supplied via either. Requires pyyaml."
    ),
)
def run_pipeline(
    genomes,
    output,
    label,
    medium,
    threads,
    diamond_bits,
    gapseq_bits,
    biomass,
    add_unique,
    soft_core,
    polish,
    engine,
    curation_db,
    pause_for_curation,
    parallel_strains,
    config_file,
):
    """Full pipeline: genomes -> pangenome -> pan-model -> strain models -> curated -> polished -> annotated -> final models."""
    import sys
    from datetime import datetime, timezone

    from parmmosahn.models import derive_strain_models
    from parmmosahn.pipeline import run_downstream, run_pan, write_pipeline_summary

    # Merge YAML config with CLI flags: CLI values win when set; for the
    # required path args, fall through to config when the CLI omitted them.
    if config_file:
        from parmmosahn.utils import load_config

        cfg = load_config(Path(config_file))
        genomes = genomes or cfg.get("genomes")
        output = output or cfg.get("output")
        label = label or cfg.get("label")
        medium = medium or cfg.get("medium")
        threads = (
            threads
            if threads != _DEFAULT_THREADS
            else cfg.get("threads", _DEFAULT_THREADS)
        )
        diamond_bits = (
            diamond_bits if diamond_bits != 150 else cfg.get("diamond_bits", 150)
        )
        gapseq_bits = gapseq_bits if gapseq_bits != 150 else cfg.get("gapseq_bits", 150)
        biomass = biomass if biomass != "pos" else cfg.get("biomass", "pos")
        engine = engine if engine != "docker" else cfg.get("engine", "docker")
        add_unique = add_unique or cfg.get("add_unique", False)
        if soft_core is None:
            soft_core = cfg.get("soft_core")
        polish = polish or cfg.get("polish", False)
        pause_for_curation = pause_for_curation or cfg.get("pause_for_curation", False)
        if parallel_strains is None:
            parallel_strains = cfg.get("parallel_strains")
        if curation_db is None:
            curation_db = cfg.get("curation_db")

    # Runtime validation — required values may come from CLI or config.
    missing = [
        name
        for name, value in (
            ("genomes", genomes),
            ("output", output),
            ("label", label),
            ("medium", medium),
        )
        if not value
    ]
    if missing:
        raise click.UsageError(
            "Missing required arguments (provide via CLI or --config): "
            + ", ".join(missing)
        )

    out_path = Path(output)
    out_path.mkdir(parents=True, exist_ok=True)

    # Set up error log file (WARNING+ only)
    error_log = out_path / "errors.log"
    error_handler = logging.FileHandler(str(error_log), mode="w")
    error_handler.setLevel(logging.WARNING)
    error_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logging.getLogger().addHandler(error_handler)

    # Set up full run log (INFO+) — append so a resumed run continues the same file
    run_log = out_path / "run.log"
    run_log_handler = logging.FileHandler(str(run_log), mode="a")
    run_log_handler.setLevel(logging.INFO)
    run_log_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logging.getLogger().addHandler(run_log_handler)

    def _teardown_log_handlers() -> None:
        """Detach and close the run/error log handlers (idempotent).

        Must run on every exit path, including the pause/abort early returns,
        or repeated in-process invocations (tests, embedding) leak open file
        handles and accumulate duplicate handlers that multiply log lines.
        """
        for _h in (run_log_handler, error_handler):
            logging.getLogger().removeHandler(_h)
            _h.close()

    t_start = time.time()
    started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    command = " ".join(sys.argv)
    steps = {}

    # Define parameters early so incremental summary writes have full context.
    # Curation DB and add_unique are included so provenance picks up the
    # actual values used (sha256 of DB, threshold actually applied).
    parameters = {
        "genomes": str(genomes),
        "medium": str(medium),
        "label": label,
        "threads": threads,
        "diamond_bits": diamond_bits,
        "gapseq_bits": gapseq_bits,
        "biomass": biomass,
        "soft_core": soft_core,
        "add_unique": add_unique,
        "engine": engine,
        "polish": polish,
        "curation_db": str(curation_db) if curation_db else None,
        "parallel_strains": parallel_strains,
        "pause_for_curation": pause_for_curation,
    }

    # --- Pre-flight validation ---
    from parmmosahn.utils import validate_input_genomes

    validate_input_genomes(Path(genomes))

    # --- Steps 1-2: Pangenome & Pan-model ---
    _echo_and_log(f"Input:  {genomes}")
    _echo_and_log(f"Output: {out_path}")
    _echo_and_log(f"Label:  {label}")
    _echo_and_log("")

    _echo_and_log("Steps 1-2: Building pangenome and pan-model...")
    t_step = time.time()
    rxn_table, fasta = run_pan(
        Path(genomes),
        out_path,
        label,
        threads=threads,
        add_unique=add_unique,
        soft_core_override=soft_core,
    )
    steps["01-pangenome + 02-panmodel"] = {
        "status": "completed",
        "elapsed": _format_elapsed(time.time() - t_step),
    }
    write_pipeline_summary(
        out_path,
        steps,
        parameters,
        command=command,
        started_at=started_at,
        total_elapsed=time.time() - t_start,
    )

    # --- Step 3: Strain models ---
    _echo_and_log("Step 3: Deriving strain-specific models...")
    t_step = time.time()
    strains_out = out_path / "03-strain_models"
    draft_models_dir = derive_strain_models(
        genomes_dir=Path(genomes),
        pangenome_fasta=fasta,
        panmodel_rxn_table=rxn_table,
        output_dir=strains_out,
        medium_file=Path(medium),
        threads=threads,
        diamond_bits=diamond_bits,
        gapseq_bits=gapseq_bits,
        biomass=biomass,
        parallel_strains=parallel_strains,
    )
    steps["03-strain_models"] = {
        "status": "completed",
        "elapsed": _format_elapsed(time.time() - t_step),
    }
    write_pipeline_summary(
        out_path,
        steps,
        parameters,
        command=command,
        started_at=started_at,
        total_elapsed=time.time() - t_start,
    )

    # --- Optional pause for human-in-the-loop curation ---
    # When stdin is a TTY, this prompts the user to edit the curation template
    # and continues inline once they hit a key. When stdin is non-interactive
    # (CI, `conda run`, nohup, tmux send-keys) the pipeline exits with a resume
    # command instead — keeping the message in `run.log` rather than only on
    # stdout, since `conda run` is known to buffer/drop stdout.
    effective_curation_db: Path | None = Path(curation_db) if curation_db else None

    if pause_for_curation:
        from parmmosahn.curation import generate_curation_template

        logger.info("Preparing curation template (running memote on drafts)...")
        template_path = generate_curation_template(draft_models_dir, out_path)
        logger.info("=== Pipeline paused for Human-in-the-Loop Curation ===")
        logger.info(
            f"ACTION REQUIRED: Fill out the curation template at {template_path.resolve()}"
        )

        if _stdin_is_tty():
            try:
                click.pause(
                    info=(
                        "Press any key once you've finished editing the "
                        "curation template to continue with Steps 4-7 "
                        "(or Ctrl-C to abort and resume later)... "
                    )
                )
            except (click.Abort, KeyboardInterrupt):
                logger.info(
                    f"Aborted. Resume with: parmmosahn project --resume -o "
                    f"{out_path.resolve()}"
                )
                _teardown_log_handlers()
                return
            logger.info("Continuing with downstream processing (Steps 4-7)...")
            effective_curation_db = template_path
        else:
            logger.info(
                f"Resume with: parmmosahn project --resume -o {out_path.resolve()}"
            )
            _teardown_log_handlers()
            return

    # --- Steps 4-7: Downstream processing ---
    logger.info(
        "Steps 4-7: Curation, polishing, annotation, gather, final evaluation..."
    )
    t_step = time.time()
    ds_result = run_downstream(
        draft_models_dir=draft_models_dir,
        output_dir=out_path,
        curation_db=effective_curation_db,
        no_polish=not polish,
        engine=engine,
    )
    steps["04-07 downstream"] = _downstream_step_entry(ds_result, time.time() - t_step)

    total_elapsed = time.time() - t_start

    # Final summary (includes model_statistics from 06-final_models)
    write_pipeline_summary(
        out_path,
        steps,
        parameters,
        command=command,
        started_at=started_at,
        total_elapsed=total_elapsed,
    )

    # Remove empty errors.log; always keep run.log
    _teardown_log_handlers()
    if error_log.exists() and error_log.stat().st_size == 0:
        error_log.unlink()

    _echo_and_log(f"\nPipeline complete! ({_format_elapsed(total_elapsed)} total)")
    _echo_and_log(f"Results: {out_path}")


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
    default=_DEFAULT_THREADS,
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
    "--parallel-strains",
    default=None,
    type=int,
    help="Number of strains to process concurrently. Default: threads//4 (auto). Each worker peaks at ~1.5 GB RAM.",
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
    parallel_strains,
):
    """Step 3: Derive strain-specific models using diamond and gapseq."""
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
        diamond_bits=diamond_bits,
        gapseq_bits=gapseq_bits,
        biomass=biomass,
        parallel_strains=parallel_strains,
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
    help="Directory containing the curated SBML models.",
)
@click.option(
    "-o",
    "--output",
    default="05b-polished_models",
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
    """Step 5b: Polish models using Palsson lab ModelPolisher (requires a container engine)."""
    from parmmosahn.curation import run_model_polisher

    models_path = Path(models)
    xml_count = len(list(models_path.rglob("*.xml")))
    click.echo(f"Input:  {models_path} ({xml_count} models)")
    click.echo(f"Output: {output}")
    click.echo(f"Engine: {engine}")

    t_start = time.time()
    run_model_polisher(models, output, engine=engine)
    click.echo(f"Done ({_format_elapsed(time.time() - t_start)} total)")


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
    default=_DEFAULT_THREADS,
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
    "--polish",
    is_flag=True,
    default=False,
    help="Enable ModelPolisher (requires Docker/Podman/Apptainer). Off by default.",
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
    "--pause-for-curation",
    "--pause_for_curation",
    "pause_for_curation",
    is_flag=True,
    help="Pause the pipeline before automatic curation to allow manual evaluation/edits.",
)
@click.option(
    "--parallel-strains",
    default=None,
    type=int,
    help="Number of strains to process concurrently. Default: threads//4 (auto). Each worker peaks at ~1.5 GB RAM.",
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
    help="Optional path to draft models for zooming/resuming. Defaults to <output>/03-strain_models/models.",
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
    polish,
    engine,
    pause_for_curation,
    parallel_strains,
    resume,
    drafts,
    config_file,
):
    """
    Rapidly project new strains (.faa or .gbff format) onto an existing pan-model and continue with downstream curation.
    This bypasses Steps 1 & 2 (PPanGGOLiN).
    Can also resume the project using the --resume flag.
    """
    import sys
    from datetime import datetime, timezone

    from parmmosahn.curation import generate_curation_template
    from parmmosahn.models import derive_strain_models
    from parmmosahn.pipeline import run_downstream, write_pipeline_summary

    if config_file:
        from parmmosahn.utils import load_config

        cfg = load_config(Path(config_file))
        genomes = genomes or cfg.get("genomes")
        fasta = fasta or cfg.get("fasta")
        panmodel_rxns = panmodel_rxns or cfg.get("panmodel_rxns")
        medium = medium or cfg.get("medium")
        threads = (
            threads
            if threads != _DEFAULT_THREADS
            else cfg.get("threads", _DEFAULT_THREADS)
        )
        diamond_bits = (
            diamond_bits if diamond_bits != 150 else cfg.get("diamond_bits", 150)
        )
        gapseq_bits = gapseq_bits if gapseq_bits != 150 else cfg.get("gapseq_bits", 150)
        biomass = biomass if biomass != "pos" else cfg.get("biomass", "pos")
        engine = engine if engine != "docker" else cfg.get("engine", "docker")
        # Downstream / Step-3 knobs newly supported in the v0.3 init-config
        # template — only override the CLI default if the user didn't pass them.
        polish = polish or cfg.get("polish", False)
        pause_for_curation = pause_for_curation or cfg.get("pause_for_curation", False)
        if parallel_strains is None:
            parallel_strains = cfg.get("parallel_strains")
        if curation_db is None:
            curation_db = cfg.get("curation_db")

    out_path = Path(output)
    out_path.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    command = " ".join(sys.argv)

    # Full run log — append so this integrates with any prior parmmosahn run log
    run_log = out_path / "run.log"
    run_log_handler = logging.FileHandler(str(run_log), mode="a")
    run_log_handler.setLevel(logging.INFO)
    run_log_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logging.getLogger().addHandler(run_log_handler)

    # Handle Resume mode
    if resume:
        # --resume runs only Steps 4-7, so any Step-3 inputs the user passed are
        # silently inert — warn rather than letting a mistaken invocation look
        # like it took effect.
        ignored = [
            flag
            for flag, val in (
                ("--genomes", genomes),
                ("--fasta", fasta),
                ("--panmodel-rxns", panmodel_rxns),
                ("--medium", medium),
                ("--parallel-strains", parallel_strains),
            )
            if val is not None
        ]
        if ignored:
            _echo_and_log(
                "Warning: --resume runs only Steps 4-7; ignoring Step-3 "
                f"option(s): {', '.join(ignored)}"
            )

        if not drafts:
            drafts = out_path / "03-strain_models" / "models"
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

        _echo_and_log(
            f"Resuming pipeline using defined curation database: {curation_db.name}"
        )

        t_step = time.time()
        ds_result = run_downstream(
            draft_models_dir=drafts,
            output_dir=out_path,
            curation_db=curation_db,
            no_polish=not polish,
            engine=engine,
        )
        write_pipeline_summary(
            out_path,
            steps={
                "04-07 downstream (resumed)": _downstream_step_entry(
                    ds_result, time.time() - t_step
                )
            },
            parameters={
                "output": str(out_path),
                "curation_db": str(curation_db),
                "medium": str(medium) if medium else None,
                "polish": polish,
                "engine": engine,
            },
            command=command,
            started_at=started_at,
            total_elapsed=time.time() - t_start,
        )
        logging.getLogger().removeHandler(run_log_handler)
        run_log_handler.close()
        return

    # Normal mode (Not resuming)
    if not all([genomes, fasta, panmodel_rxns, medium]):
        raise click.UsageError(
            "Missing required arguments for normal project strain derivation."
        )

    # 1. Run Strain Derivation (Step 3)
    click.echo(step_banner(3))
    strains_out = out_path / "03-strain_models"

    # Defensive runtime validation for required CLI args; do not use assert here,
    # because assertions can be removed with optimized Python execution.
    if genomes is None or fasta is None or panmodel_rxns is None or medium is None:
        raise click.UsageError(
            "Missing required arguments: genomes, fasta, panmodel_rxns, and medium must be provided"
        )

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
        parallel_strains=parallel_strains,
    )

    # See run_pipeline for the rationale on logger.info over click.echo and the
    # TTY-aware inline-continue path (B7 in docs/ROADMAP.md).
    effective_curation_db: Path | None = Path(curation_db) if curation_db else None

    if pause_for_curation:
        logger.info("=== Preparing Human-in-the-Loop Curation (running memote) ===")
        template_path = generate_curation_template(draft_models_dir, out_path)
        logger.info("=== Pipeline paused for Human-in-the-Loop Curation ===")
        logger.info(
            f"ACTION REQUIRED: Please fill out the curation template at "
            f"{template_path.resolve()}"
        )

        if _stdin_is_tty():
            try:
                click.pause(
                    info=(
                        "Press any key once you've finished editing the "
                        "curation template to continue with Steps 4-7 "
                        "(or Ctrl-C to abort and resume later)... "
                    )
                )
            except (click.Abort, KeyboardInterrupt):
                logger.info(
                    f"Aborted. Resume with: parmmosahn project --resume -o "
                    f"{out_path.resolve()}"
                )
                logging.getLogger().removeHandler(run_log_handler)
                run_log_handler.close()
                return
            logger.info("Continuing with downstream processing (Steps 4-7)...")
            effective_curation_db = template_path
        else:
            logger.info(
                f"Resume with: parmmosahn project --resume -o {out_path.resolve()}"
            )
            logging.getLogger().removeHandler(run_log_handler)
            run_log_handler.close()
            return

    # 2. Pipe directly into downstream processing (Steps 4-7)
    t_downstream = time.time()
    ds_result = run_downstream(
        draft_models_dir=draft_models_dir,
        output_dir=out_path,
        curation_db=effective_curation_db,
        no_polish=not polish,
        engine=engine,
    )
    write_pipeline_summary(
        out_path,
        steps={
            "03-strain_models": {
                "status": "completed",
                "elapsed": _format_elapsed(t_downstream - t_start),
            },
            "04-07 downstream": _downstream_step_entry(
                ds_result, time.time() - t_downstream
            ),
        },
        parameters={
            "genomes": str(genomes),
            "fasta": str(fasta),
            "panmodel_rxns": str(panmodel_rxns),
            "medium": str(medium),
            "threads": threads,
            "diamond_bits": diamond_bits,
            "gapseq_bits": gapseq_bits,
            "biomass": biomass,
            "polish": polish,
            "engine": engine,
            "curation_db": str(curation_db) if curation_db else None,
            "parallel_strains": parallel_strains,
        },
        command=command,
        started_at=started_at,
        total_elapsed=time.time() - t_start,
    )
    logging.getLogger().removeHandler(run_log_handler)
    run_log_handler.close()


_INIT_CONFIG_TEMPLATE = """\
# PaRMMoSaHN configuration file
# Generated by: parmmosahn init-config
# All values shown are defaults. Remove or comment out any you don't need.
# CLI flags always override values set here.
# You can also use environment variables: PARMMOSAHN_THREADS=32, etc.

# --- Input paths (pan / run / project) ---
# genomes: /path/to/genomes          # Directory of genome files (.gbff or .faa)
# fasta: /path/to/softcore.faa       # Existing soft-core FASTA (project command)
# panmodel_rxns: /path/to/rxns.tbl   # Existing pan-model reaction table (project command)
# medium: /path/to/medium.csv        # Media CSV file for gapfilling

# --- Output (pan / run / project) ---
# output: .                          # Base output directory (creates 01-pangenome/, etc.)

# --- Pangenome (pan / run commands) ---
# label: myproject                   # Label for the pangenome
# add_unique: false                  # Include unique/cloud genes in the pan-model
# soft_core: null                    # Override soft-core threshold (fraction 0-1).
#                                    # Default: 2/n_genomes (biologically motivated).

# --- Performance (run / project / strains) ---
threads: 14
diamond_bits: 150
gapseq_bits: 150
biomass: pos                         # Biomass reaction type: 'pos' or 'neg'
# parallel_strains: null             # Strains processed in parallel in Step 3.
#                                    # null = auto (threads // 4); each worker
#                                    # peaks at ~1.5 GB RAM. Lower this if
#                                    # memote OOMs under WSL2.

# --- Downstream (run / project) ---
polish: false                        # Run ModelPolisher (requires Docker/Podman/Apptainer).
                                     # Off by default — the bundled v2.1-beta
                                     # has a fragile MIRIAM/identifiers.org
                                     # network fetch at startup.
engine: docker                       # Container engine when polish: true.
# curation_db: /path/to/curation.xlsx  # Existing curation Excel. If set, Step
#                                      # 4a (memote on drafts → template
#                                      # generation) is skipped and Step 4b
#                                      # applies this DB directly.
# pause_for_curation: false          # If true, Steps 1-3 + 4a run, then the
#                                    # pipeline pauses for manual curation
#                                    # edits before continuing to 4b-7.
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
def annotate(models, output, rename_to_bigg):
    """Enrich model annotations for FAIR compliance (genes, reactions)."""
    from parmmosahn.annotation import enrich_models_directory

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
