import contextlib
import functools
import json
import logging
import multiprocessing
import os
import re
import subprocess
import time
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any, Dict, Tuple, cast

import cobra as cb
import polars as pl
from cobra.core import Model, Reaction

# Import memote API directly
from memote.suite import api as memote_api
from memote.suite.reporting import ReportConfiguration

from parmmosahn.constants import step_banner
from parmmosahn.utils import ToolExecutionError

logger = logging.getLogger(__name__)

# ==========================================
# Phase 1: memote evaluation & DB Generation
# ==========================================


def _run_memote_snapshot_inprocess(
    model_path: Path, output_json: Path, output_html: Path
) -> None:
    """Run memote snapshot in the current process (high memory usage)."""
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    # cobra's SBML reader logs one INFO/WARNING line per boundary metabolite
    # ("Adding exchange reaction EX_..."), which floods the J9 log with
    # hundreds of lines per model and buries the per-strain progress. These
    # are informational, not actionable. Silence the cobra logger here (this
    # runs in the spawned memote child, so it does not affect the parent's
    # logging) before validate_model triggers the SBML read.
    for _name in ("cobra", "cobra.io.sbml", "cobra.core"):
        logging.getLogger(_name).setLevel(logging.ERROR)

    model_obj, sbml_ver, notifications = memote_api.validate_model(str(model_path))
    if model_obj is None:
        logger.error(
            f"Failed to load model {model_path.name} for memote due to SBML errors."
        )
        # Raise rather than return: in the spawned child this yields a non-zero
        # exit code, and in the in-process path it is caught by the caller. Both
        # then record the model as FAILED instead of silently counting it as a
        # passed snapshot and dropping it from the curation template / report.
        raise ToolExecutionError(
            f"memote could not load {model_path.name} (SBML errors)."
        )

    pytest_args = ["--tb", "no", "-q"]

    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            output = memote_api.test_model(
                model_obj, sbml_version=sbml_ver, results=True, pytest_args=pytest_args
            )
            code, results = cast(Tuple[int, Any], output)

    config = ReportConfiguration.load()

    html_report = memote_api.snapshot_report(results, config=config, html=True)
    if html_report:
        # Write atomically: a kill mid-write must never leave a half-formed file
        # that downstream consumers (or our own --skip-existing logic) treat as valid.
        tmp_html = output_html.with_suffix(output_html.suffix + ".tmp")
        with tmp_html.open("w", encoding="utf-8") as f:
            f.write(cast(str, html_report))
        tmp_html.replace(output_html)

    json_report = memote_api.snapshot_report(results, config=config, html=False)
    if json_report:
        tmp_json = output_json.with_suffix(output_json.suffix + ".tmp")
        with tmp_json.open("w", encoding="utf-8") as f:
            f.write(cast(str, json_report))
        tmp_json.replace(output_json)


def run_memote_snapshot(
    model_path: Path,
    output_json: Path,
    output_html: Path,
    subprocess_mode: bool = True,
    timeout: int = 900,
) -> bool:
    """
    Runs memote report snapshot on a single model.

    Parameters
    ----------
    subprocess_mode : bool
        If True (default), run memote in a separate process so that memory
        is fully reclaimed between models. This prevents OOM kills when
        evaluating multiple large models sequentially.
    timeout : int
        Maximum seconds to wait for memote per model (default: 900 = 15 min).
    """
    logger.debug(f"Running memote snapshot on {model_path.name}...")

    if not subprocess_mode:
        try:
            _run_memote_snapshot_inprocess(model_path, output_json, output_html)
        except Exception as exc:
            logger.error(f"Memote snapshot failed for {model_path.name}: {exc}")
            _cleanup_partial_outputs(output_json, output_html)
            return False
        logger.debug(f"Memote snapshot completed for {model_path.name}.")
        return True

    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    # A plain ``spawn`` Process, not ProcessPoolExecutor. Two deadlocks were
    # observed when the executor variant ran inside the ThreadPoolExecutor that
    # ``parmmosahn evaluate --workers N`` uses (deadlock_capture_20260517):
    #
    #   1. ProcessPoolExecutor defaults to the ``fork`` start method. Forking a
    #      multi-threaded process (the outer thread pool + jemalloc bg threads)
    #      wedged the child at bootstrap in multiprocessing/util.py
    #      _flush_std_streams, holding a stdio lock copied from the parent.
    #   2. Even when future.result(timeout=) raised, leaving the
    #      ``with ProcessPoolExecutor`` block called shutdown(wait=True), which
    #      os.waitpid()-blocked forever on the still-running hung memote child.
    #      The timeout cancelled the *future*, never the *process*.
    #
    # ``spawn`` gives the child a fresh interpreter with no inherited locks
    # (kills cause 1). join(timeout) + terminate() + kill() guarantees the
    # child is reaped within ``timeout`` no matter how memote hangs (kills
    # cause 2): a pathological model can never freeze the run.
    ctx = multiprocessing.get_context("spawn")
    proc = ctx.Process(
        target=_run_memote_snapshot_inprocess,
        args=(model_path, output_json, output_html),
        name=f"memote-{model_path.stem}",
    )
    proc.start()
    proc.join(timeout)

    if proc.is_alive():
        logger.error(
            f"Memote exceeded {timeout}s for {model_path.name}; killing and "
            "skipping (run continues, --skip-existing will retry)."
        )
        proc.terminate()
        proc.join(30)
        if proc.is_alive():
            proc.kill()
            proc.join()
        _cleanup_partial_outputs(output_json, output_html)
        return False

    if proc.exitcode != 0:
        logger.error(
            f"Memote subprocess for {model_path.name} exited with "
            f"{proc.exitcode}; skipping."
        )
        _cleanup_partial_outputs(output_json, output_html)
        return False

    logger.debug(f"Memote snapshot completed for {model_path.name}.")
    return True


def _cleanup_partial_outputs(output_json: Path, output_html: Path) -> None:
    """Remove half-written ``.tmp`` files left by a killed memote child.

    ``_run_memote_snapshot_inprocess`` writes atomically (``.tmp`` then
    ``replace``), so a kill never produces a corrupt final file; it can only
    leave an orphan ``.tmp``. Removing it keeps the output dir clean so
    ``--skip-existing`` sees the model as not-yet-done and retries it.
    """
    for path in (output_json, output_html):
        tmp = path.with_suffix(path.suffix + ".tmp")
        with contextlib.suppress(OSError):
            tmp.unlink(missing_ok=True)


def _excel_safe_df(df: "pl.DataFrame") -> "pl.DataFrame":
    """Prefix any string cell beginning with ``= + - @`` with a quote.

    Spreadsheet apps interpret such a leading character as a formula when the
    curation template is opened; the leading quote forces literal text (a
    CSV/XLSX-injection guard). Non-string columns are left untouched.
    """
    leading = ("=", "+", "-", "@")
    str_cols = [c for c in df.columns if df.schema[c] == pl.Utf8]
    if not str_cols:
        return df
    return df.with_columns(
        [
            pl.when(pl.col(c).str.slice(0, 1).is_in(leading))
            .then(pl.concat_str([pl.lit("'"), pl.col(c)]))
            .otherwise(pl.col(c))
            .alias(c)
            for c in str_cols
        ]
    )


def extract_memote_problems(
    report_dir: Path, model_dir: Path, output_excel: Path
) -> None:
    """
    Parses memote JSON reports to find imbalances and duplicates.
    Writes a 3-sheet Excel template matching the PaRMMoSaHN curation DB format.
    """
    logger.info("Extracting problem balances and duplicates from memote reports...")
    output_excel.parent.mkdir(parents=True, exist_ok=True)

    def read_test_results(test_key: str) -> dict:
        rxns: set[Any] = set()
        first_encountered: dict[Any, str] = {}
        for report_file in report_dir.glob("*.json"):
            model_name = report_file.stem
            try:
                with report_file.open("r", encoding="utf-8") as f:
                    cont = json.load(f)
                encountered_rxns = (
                    cont.get("tests", {}).get(test_key, {}).get("data", [])
                )
                if not encountered_rxns:
                    continue
                if isinstance(encountered_rxns[0], list):
                    encountered_rxns = [tuple(i) for i in encountered_rxns]
                new_rxns = set(encountered_rxns).difference(rxns)
                for r in new_rxns:
                    first_encountered[r] = model_name
                rxns.update(new_rxns)
            except Exception as exc:
                # Don't let one malformed/corrupt report silently empty the
                # template — a systematic issue (e.g. a memote schema change)
                # must be visible in the log.
                logger.warning(
                    f"Could not parse memote report {report_file.name} "
                    f"for {test_key}: {exc}"
                )
        return first_encountered

    first_encounter_charge = read_test_results("test_reaction_charge_balance")
    first_encounter_mass = read_test_results("test_reaction_mass_balance")
    first_encounter_dupl = read_test_results("test_find_duplicate_reactions")

    # Cache loaded models to avoid re-reading the same XML for each of its
    # problem rows, but BOUND the cache (LRU) so a 1000-model cohort cannot pin
    # hundreds of full cobra models in memory at once. maxsize trades a few
    # re-reads for a hard memory ceiling.
    @functools.lru_cache(maxsize=16)
    def get_cached_model(model_name: str):
        # Check the new subfolder architecture, then fall back to a flat layout.
        model_path = model_dir / model_name / f"{model_name}.xml"
        if not model_path.exists():
            model_path = model_dir / f"{model_name}.xml"

        if model_path.exists():
            return cast(Model, cb.io.read_sbml_model(str(model_path)))
        return None

    def get_duplicates_df() -> pl.DataFrame:
        df_rows = []
        # Anchor on the SEED reaction id (``rxn`` + its full numeric block) rather
        # than the first 5-digit run anywhere in the id: the latter mis-parses
        # 6+ digit ids and can match a ``cpd`` number inside an exchange id.
        id_regex = re.compile(r"rxn0*(\d+)")
        for rxn_pair, model_name in first_encounter_dupl.items():
            network = get_cached_model(model_name)
            if not network:
                continue
            try:
                rxn_pair_obj = [
                    cast(Reaction, network.reactions.get_by_id(str(r)))
                    for r in rxn_pair
                ]

                # No 99999 fallback: warn and sort unparseable ids to the back.
                rxn_ids: list[float] = []
                for r in rxn_pair_obj:
                    match = id_regex.search(r.id)
                    if match:
                        rxn_ids.append(float(match.group(1)))
                    else:
                        logger.warning(
                            f"Reaction ID '{r.id}' is not a SEED rxn id; "
                            f"cannot order it for duplicate sorting."
                        )
                        rxn_ids.append(
                            float("inf")
                        )  # Sort to the back safely without faking an ID

                min_idx = rxn_ids.index(min(rxn_ids))
                r1 = cast(Reaction, rxn_pair_obj.pop(min_idx))
                r2 = cast(Reaction, rxn_pair_obj[0])

                df_rows.append(
                    {
                        "Pair": str(rxn_pair),
                        "Decision": "",
                        "ID 1": r1.id,
                        "Name 1": r1.name,
                        "Equation 1": r1.build_reaction_string(),
                        "Lower bound 1": r1.lower_bound,
                        "Upper bound 1": r1.upper_bound,
                        "ID 2": r2.id,
                        "Name 2": r2.name,
                        "Equation 2": r2.build_reaction_string(),
                        "Lower bound 2": r2.lower_bound,
                        "Upper bound 2": r2.upper_bound,
                    }
                )
            except Exception as e:
                logger.debug(f"Skipping duplicate pair {rxn_pair}: {e}")

        return (
            pl.DataFrame(df_rows)
            if df_rows
            else pl.DataFrame(
                schema=["Pair", "Decision", "ID 1", "Name 1", "Equation 1"]
            )
        )

    def get_imbalance_df() -> pl.DataFrame:
        df_rows = []
        combined = [(r, m, "charge") for r, m in first_encounter_charge.items()] + [
            (r, m, "mass") for r, m in first_encounter_mass.items()
        ]

        for rxn_id, model_name, bal_type in combined:
            network = get_cached_model(model_name)
            if not network:
                continue
            try:
                rxn_obj = cast(Reaction, network.reactions.get_by_id(str(rxn_id)))
                df_rows.append(
                    {
                        "reaction_id": rxn_obj.id,
                        "name": rxn_obj.name,
                        "new_formula": rxn_obj.build_reaction_string(),
                        "balance_type": bal_type,
                        "modification": "",
                        "Comment": str(rxn_obj.check_mass_balance()).strip("{}"),
                    }
                )
            except Exception as e:
                logger.debug(f"Skipping imbalanced reaction {rxn_id}: {e}")
        return pl.DataFrame(df_rows)

    # Write DataFrames to Excel
    logger.info(f"Writing curation template to {output_excel.name}...")

    df_dupl = get_duplicates_df()
    if df_dupl.is_empty():
        df_dupl = pl.DataFrame(
            schema=["Pair", "Decision", "ID 1", "Name 1", "Equation 1"]
        )

    df_imb = get_imbalance_df()
    if df_imb.is_empty():
        df_imb = pl.DataFrame(
            schema=[
                "reaction_id",
                "name",
                "new_formula",
                "balance_type",
                "modification",
                "Comment",
            ]
        )

    df_ign = pl.DataFrame(
        schema=["reaction_id", "name", "Imbalance", "Reason to leave intact", "Comment"]
    )

    # Using xlsxwriter engine in Polars to write multiple sheets efficiently.
    # The Workbook must be created explicitly; passing a path to write_excel
    # can auto-close the workbook after the first sheet.
    from xlsxwriter import Workbook

    with Workbook(str(output_excel)) as wb:
        _excel_safe_df(df_dupl).write_excel(
            workbook=wb, worksheet="duplicate_reactions"
        )
        _excel_safe_df(df_imb).write_excel(workbook=wb, worksheet="curated_imbalances")
        df_ign.write_excel(workbook=wb, worksheet="ignored_imbalances")


# ==========================================
# Phase 2: Applying the Curation Database
# ==========================================


def _apply_duplicate_decisions(
    model: "Model",
    df_dupl: "pl.DataFrame",
    stats: dict | None = None,
) -> int:
    """Apply duplicate reaction removal decisions from the curation DB.

    Decision logic:
      0: Drop both reactions
      1: Keep ID 1, drop ID 2
      2: Keep ID 2, drop ID 1
      3: Keep both (do nothing)
      blank: Skip and continue

    When *stats* is given, each row's hit count is accumulated under the key
    ``duplicate_reactions`` keyed by ``"<ID 1>|<ID 2>"`` so the orchestrator
    can emit a per-row applicability report across strains.

    Returns:
        Number of reactions removed.
    """
    sheet_stats = None
    if stats is not None:
        sheet_stats = stats.setdefault("duplicate_reactions", {})

    changes_made = 0
    for row in df_dupl.iter_rows(named=True):
        decision = row.get("Decision")
        id1 = row.get("ID 1")
        id2 = row.get("ID 2")
        row_key = f"{id1}|{id2}"
        row_entry = (
            sheet_stats.setdefault(
                row_key,
                {
                    "id1": id1,
                    "id2": id2,
                    "decision": str(decision).strip() if decision is not None else "",
                    "strains_hit": 0,
                    "reactions_removed": 0,
                },
            )
            if sheet_stats is not None
            else None
        )

        if decision is None or str(decision).strip() == "":
            continue

        decision_str = str(decision).strip()
        to_drop: list = []

        if decision_str == "0":
            to_drop = [id1, id2]
        elif decision_str == "1":
            to_drop = [id2]
        elif decision_str == "2":
            to_drop = [id1]
        elif decision_str == "3":
            continue
        else:
            logger.warning(
                f"Unknown decision '{decision_str}' for pair {id1}/{id2}. Skipping."
            )
            continue

        hit_this_strain = False
        for rxn_id in to_drop:
            if rxn_id and rxn_id in model.reactions:
                rxn_to_drop = cast(Reaction, model.reactions.get_by_id(rxn_id))
                rxn_to_drop.remove_from_model(remove_orphans=True)
                changes_made += 1
                if row_entry is not None:
                    row_entry["reactions_removed"] += 1
                hit_this_strain = True
        if row_entry is not None and hit_this_strain:
            row_entry["strains_hit"] += 1

    return changes_made


def _apply_formula_fixes(
    model: "Model",
    df_cur: "pl.DataFrame",
    model_name: str,
    stats: dict | None = None,
) -> int:
    """Apply stoichiometric formula corrections from the curated_imbalances sheet.

    When *stats* is given, each row's hit count is accumulated under
    ``curated_imbalances`` keyed by ``reaction_id`` so the orchestrator can
    surface DB rows that fired on zero strains (i.e. stale corrections).

    Returns:
        Number of reactions updated.
    """
    sheet_stats = None
    if stats is not None:
        sheet_stats = stats.setdefault("curated_imbalances", {})

    changes_made = 0
    for row in df_cur.iter_rows(named=True):
        rxn_id = row.get("reaction_id")
        new_formula = row.get("new_formula")
        row_entry = (
            sheet_stats.setdefault(
                str(rxn_id) if rxn_id is not None else "",
                {
                    "reaction_id": rxn_id,
                    "new_formula": new_formula,
                    "strains_hit": 0,
                },
            )
            if sheet_stats is not None and rxn_id is not None
            else None
        )

        if rxn_id is None or new_formula is None or str(rxn_id) not in model.reactions:
            continue

        rxn = cast(Reaction, model.reactions.get_by_id(str(rxn_id)))
        old_formula = rxn.build_reaction_string()

        if str(new_formula).strip() != old_formula:
            try:
                rxn.build_reaction_from_string(str(new_formula).strip())
                changes_made += 1
                if row_entry is not None:
                    row_entry["strains_hit"] += 1
            except Exception as e:
                logger.warning(
                    f"Could not apply new formula for {rxn_id} in {model_name}: {e}"
                )

    return changes_made


def apply_curations(
    model_path: Path,
    curation_db_path: Path,
    output_dir: Path,
    stats: dict | None = None,
) -> Path | None:
    """Read the curation database and apply all fixes to a specific strain model.

    Orchestrates duplicate-reaction decisions and formula corrections. The public
    signature is unchanged for callers that don't care about per-row stats;
    pass *stats* (a shared dict) to accumulate per-row hit counts across
    strains for downstream applicability reporting.

    Skip-existing: if the output already exists it is returned untouched so a
    resumed run does not re-curate. The output is written atomically, so an
    existing file is always complete (never a truncated partial). To re-apply an
    updated curation database, clear the output directory first.
    """
    out_file = output_dir / model_path.name

    if out_file.exists():
        return out_file

    try:
        model = cast(Model, cb.io.read_sbml_model(str(model_path)))
    except Exception as e:
        logger.error(f"Failed to load {model_path.name} for curation: {e}")
        return None

    try:
        df_dupl = pl.read_excel(
            curation_db_path, sheet_name="duplicate_reactions", engine="openpyxl"
        )
        # Accept both the canonical sheet name and the legacy name for balance
        # corrections. "curated_imbalances" is preferred; "corrected_balances"
        # is the fallback for older curation databases.
        try:
            df_cur = pl.read_excel(
                curation_db_path, sheet_name="curated_imbalances", engine="openpyxl"
            )
        except Exception:
            df_cur = pl.read_excel(
                curation_db_path, sheet_name="corrected_balances", engine="openpyxl"
            )
    except Exception as e:
        logger.error(f"Could not read curation database {curation_db_path}: {e}")
        return None

    changes_made = 0
    changes_made += _apply_duplicate_decisions(model, df_dupl, stats=stats)
    changes_made += _apply_formula_fixes(model, df_cur, model_path.name, stats=stats)

    # Atomic write: a crash mid-write must not leave a truncated SBML that the
    # skip-existing guard above would later mistake for a completed curation.
    tmp_out = out_file.with_suffix(out_file.suffix + ".tmp")
    cb.io.write_sbml_model(model, str(tmp_out))
    os.replace(tmp_out, out_file)
    logger.debug(f"Applied {changes_made} curations to {model_path.name}")

    return out_file


def write_curation_application_report(
    stats: dict,
    n_strains: int,
    output_path: Path,
) -> Path:
    """Write a per-row applicability report from accumulated curation stats.

    Each row of the curation DB is emitted with the number of strains in which
    its rule fired. Rows with ``strains_hit == 0`` are warned about — they
    represent DB rules that have no targets in the current cohort (typically
    stale rules from an earlier gap-fill iteration; see A9 §F).

    Args:
        stats: dict produced by :func:`apply_curations` across all strains.
        n_strains: total strains processed (denominator for the report).
        output_path: target TSV path.

    Returns:
        Path to the written TSV.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[str] = [
        "sheet\tidentifier\tdetail\tstrains_hit\tn_strains\tfraction\tnotes"
    ]
    zero_hit_count = 0

    dupl = stats.get("duplicate_reactions", {})
    for key, entry in sorted(dupl.items()):
        hit = int(entry.get("strains_hit", 0))
        frac = (hit / n_strains) if n_strains else 0.0
        note = "" if hit > 0 else "no targets in cohort"
        if hit == 0:
            zero_hit_count += 1
        detail = f"decision={entry.get('decision', '')}"
        rows.append(
            f"duplicate_reactions\t{key}\t{detail}\t"
            f"{hit}\t{n_strains}\t{frac:.4f}\t{note}"
        )

    cur = stats.get("curated_imbalances", {})
    for key, entry in sorted(cur.items()):
        hit = int(entry.get("strains_hit", 0))
        frac = (hit / n_strains) if n_strains else 0.0
        note = "" if hit > 0 else "no targets in cohort"
        if hit == 0:
            zero_hit_count += 1
        nf = entry.get("new_formula") or ""
        # Keep TSV one-line-per-row; truncate runaway formulas.
        nf_short = (str(nf)[:120] + "…") if len(str(nf)) > 120 else str(nf)
        rows.append(
            f"curated_imbalances\t{key}\tnew_formula={nf_short}\t"
            f"{hit}\t{n_strains}\t{frac:.4f}\t{note}"
        )

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    n_dupl = len(dupl)
    n_cur = len(cur)
    logger.info(
        f"Curation applicability report written to {output_path} "
        f"({n_dupl} duplicate rows + {n_cur} imbalance rows across {n_strains} strains)"
    )
    if zero_hit_count:
        logger.warning(
            f"{zero_hit_count} curation DB row(s) fired on 0/{n_strains} strains — "
            f"likely stale rules from an earlier gap-fill iteration "
            f"(see {output_path.name})."
        )
    return output_path


def _await_polisher_databases(
    base_cmd: list[str],
    compose_file: Path,
    compose_dir: Path,
    run_env: dict,
    services: tuple[str, ...] = ("biggdb", "adb"),
    timeout: int = 300,
) -> None:
    """Block until the ModelPolisher PostgreSQL backends accept connections.

    Polls ``pg_isready`` per service rather than sleeping a fixed interval, so
    fast hosts proceed at once and slow ones still get a generous budget. If the
    probe can't run (engine lacks ``exec``, or no ``pg_isready`` in the image)
    the poll simply runs to *timeout* and then returns, matching the old
    fixed-wait behaviour. Factored out so it can be stubbed in tests.
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        all_ready = True
        for svc in services:
            probe = base_cmd + [
                "-f",
                str(compose_file),
                "exec",
                "-T",
                svc,
                "pg_isready",
            ]
            try:
                rc = subprocess.run(
                    probe,
                    cwd=str(compose_dir),
                    env=run_env,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    timeout=30,
                ).returncode
            except (subprocess.SubprocessError, OSError):
                rc = 1
            if rc != 0:
                all_ready = False
                break
        if all_ready:
            logger.info("ModelPolisher databases ready.")
            return
        time.sleep(5)
    logger.warning(
        f"Database readiness probe did not succeed within {timeout}s; "
        "proceeding anyway (ModelPolisher may fail if a DB is not up)."
    )


def run_model_polisher(
    model_dir: str | Path,
    output_dir: str | Path,
    engine: str = "docker",
    timeout: int = 1800,
) -> Path:
    """
    Executes the Palsson lab ModelPolisher via Docker Compose or alternative engines.
    Automatically provisions the local BiGG and Annotation databases.

    *timeout* bounds each per-model polish (default 1800 s = 30 min) so a single
    wedged container cannot stall the whole batch; a model that exceeds it is
    logged, skipped, and counted as a failure.
    """
    model_dir = Path(model_dir).resolve()
    output_dir = Path(output_dir).resolve()

    if engine == "docker":
        base_cmd = ["docker", "compose"]
    elif engine == "podman":
        base_cmd = ["podman-compose"]
    elif engine == "apptainer":
        base_cmd = ["apptainer", "compose"]
    else:
        raise ValueError(f"Unsupported container engine: {engine}")

    models_out = output_dir / "models"
    logs_out = output_dir / "logs"
    models_out.mkdir(parents=True, exist_ok=True)
    logs_out.mkdir(parents=True, exist_ok=True)

    xml_models = list(model_dir.rglob("*.xml"))
    if not xml_models:
        logger.warning(f"No XML models found in {model_dir} — skipping polish.")
        return models_out
    logger.info(f"Polishing {len(xml_models)} models with ModelPolisher ({engine})...")

    run_env = os.environ.copy()
    run_env["MODEL_DIR"] = str(model_dir)
    run_env["OUTPUT_DIR"] = str(models_out)
    run_env["LOGS_DIR"] = str(logs_out)

    compose_resource = files("parmmosahn.data").joinpath("docker-compose.yml")

    with as_file(compose_resource) as compose_file:
        compose_dir = compose_file.parent

        # STEP 1: Boot up the background databases
        logger.info(
            "Starting local ModelPolisher databases (this may take a moment)..."
        )
        subprocess.run(
            base_cmd
            + ["-f", str(compose_file), "up", "-d", "biggdb", "adb", "validator"],
            cwd=str(compose_dir),
            env=run_env,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Wait until the PostgreSQL backends actually accept connections —
        # ModelPolisher fails if it connects too early.
        logger.info("Waiting for ModelPolisher databases to accept connections...")
        _await_polisher_databases(base_cmd, compose_file, compose_dir, run_env)

        # Downgrades SBML headers from L3V2 to L3V1 for compatibility with the beta version of ModelPolisher. This is a temporary workaround until we can upgrade to the latest ModelPolisher release that supports L3V2.
        logger.info("Downgrading SBML headers to L3V1 for beta compatibility...")
        for xml in xml_models:
            content = xml.read_text(encoding="utf-8")
            if 'level="3" version="2"' in content:
                content = content.replace(
                    'level="3" version="2"', 'level="3" version="1"'
                )
                content = content.replace(
                    "level3/version2/core", "level3/version1/core"
                )
                xml.write_text(content, encoding="utf-8")
        # -------------------------------------------------------

        success_count = 0
        fail_count = 0

        try:
            # STEP 2: Process the models
            n_models = len(xml_models)
            for i, xml in enumerate(xml_models, 1):
                strain_name = xml.stem
                logger.info(f"[{i}/{n_models}] Polishing {strain_name}...")

                cmd = base_cmd + [
                    "-f",
                    str(compose_file),
                    "run",
                    "-T",
                    "--rm",
                    "polisher",
                    "java",
                    "-jar",
                    "ModelPolisher-2.1-beta.jar",
                    f"--input=/models/{xml.name}",
                    f"--output=/output/{xml.name}",
                    "--annotate-with-bigg=true",
                ]

                # Stream child output straight to the per-model log file rather
                # than buffering it in RAM — a verbose/looping polish across
                # 1000 models would otherwise accumulate unbounded output.
                log_file = logs_out / f"{strain_name}.log"
                try:
                    with log_file.open("w", encoding="utf-8") as log_fh:
                        subprocess.run(
                            cmd,
                            cwd=str(compose_dir),
                            env=run_env,
                            stdout=log_fh,
                            stderr=subprocess.STDOUT,
                            text=True,
                            check=True,
                            timeout=timeout,
                        )
                    success_count += 1

                except subprocess.TimeoutExpired:
                    logger.error(
                        f"ModelPolisher timed out after {timeout}s on "
                        f"{strain_name}; skipping (see {log_file})."
                    )
                    fail_count += 1
                    continue

                except subprocess.CalledProcessError:
                    logger.error(
                        f"ModelPolisher failed on {strain_name} "
                        f"(see {log_file} for details)."
                    )
                    fail_count += 1
                    continue

        finally:
            # STEP 3: Shut down databases
            logger.info("Shutting down local ModelPolisher databases...")
            subprocess.run(
                base_cmd + ["-f", str(compose_file), "down"],
                cwd=str(compose_dir),
                env=run_env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    if fail_count > 0:
        logger.warning(
            f"Polishing finished with {fail_count} failures and {success_count} successes."
        )
    else:
        logger.info(f"Successfully polished all {success_count} models.")

    return models_out


# ==========================================
# Phase 3: Hidden Evaluators (Pipeline Auto-Triggers)
# ==========================================


def generate_curation_template(model_dir: Path, output_dir: Path) -> Path:
    """Hidden Evaluator: Runs memote on drafts and generates the Excel template."""
    reports_dir = output_dir / "07-memote_reports" / "draft"
    reports_dir.mkdir(parents=True, exist_ok=True)

    models = sorted(model_dir.rglob("*.xml"))
    n = len(models)
    logger.info(f"Running automated memote evaluation on {n} draft models...")

    failed = []
    for i, model in enumerate(models, 1):
        logger.info(f"  [{i}/{n}] Evaluating {model.stem}...")
        out_json = reports_dir / f"{model.stem}.json"
        out_html = reports_dir / f"{model.stem}.html"
        if not run_memote_snapshot(model, out_json, out_html):
            failed.append(model.stem)

    if failed:
        logger.warning(f"Memote failed for {len(failed)}/{n} models: {failed}")

    out_excel = output_dir / "curation_template.xlsx"
    extract_memote_problems(reports_dir, model_dir, out_excel)
    return out_excel


def evaluate_final_models(
    model_dir: Path, output_dir: Path, skip_existing: bool = True
) -> None:
    """Hidden Evaluator: Runs memote on final models to prove FAIR compliance.

    With *skip_existing* (default), a model whose memote JSON report already
    exists and is non-empty is skipped, so a resumed run does not repeat the
    slow per-model evaluation. A fresh run finds an empty reports dir and so
    evaluates everything.
    """
    reports_dir = output_dir / "07-memote_reports" / "final"
    reports_dir.mkdir(parents=True, exist_ok=True)

    models = sorted(model_dir.rglob("*.xml"))
    if not models:
        logger.warning("No final XML models found to evaluate.")
        return

    n = len(models)
    logger.info(step_banner(7))
    logger.info(f"Running final memote evaluation on {n} polished models...")
    failed = []
    skipped = 0
    for i, model in enumerate(models, 1):
        out_json = reports_dir / f"{model.stem}.json"
        out_html = reports_dir / f"{model.stem}.html"
        if skip_existing and out_json.exists() and out_json.stat().st_size > 0:
            skipped += 1
            continue
        logger.info(f"  [{i}/{n}] Evaluating {model.stem}...")
        if not run_memote_snapshot(model, out_json, out_html):
            failed.append(model.stem)

    if skipped:
        logger.info(f"Skipped {skipped}/{n} model(s) with existing reports.")
    if failed:
        logger.warning(f"Memote failed for {len(failed)}/{n} models: {failed}")
    logger.info(f"Final FAIR-compliance reports saved to {reports_dir}")
