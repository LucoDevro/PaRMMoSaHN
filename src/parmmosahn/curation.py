import contextlib
import functools
import json
import logging
import multiprocessing
import os
import re
import warnings
from pathlib import Path
from typing import Any, Tuple, cast

import cobra as cb
import polars as pl
from cobra.core import Model, Reaction

# Import memote API directly
from memote.suite import api as memote_api
from memote.suite.reporting import ReportConfiguration

from parmmosahn.utils import ToolExecutionError

logger = logging.getLogger(__name__)


# ==========================================
# Phase 0: create an empty curation DB
# ==========================================


def init_curation_table(output_excel: Path, force: bool = False) -> None:
    """
    Write an empty PaRMMoSaHN curation DB.
    """
    output_excel.parent.mkdir(parents=True, exist_ok=True)
    if output_excel.is_file() and not force:
        raise ToolExecutionError(
            'Output file already exists! Use --force to overwrite it.'
            )
    
    df_dupl = pl.DataFrame(
        schema=["Pair", "Decision",
                "ID 1", "Name 1", "Equation 1", "Lower bound 1", "Upper bound 1",
                "ID 2", "Name 2", "Equation 2", "Lower bound 2", "Upper bound 2"]
        )
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
    
    from xlsxwriter import Workbook

    with Workbook(str(output_excel)) as wb:
        _excel_safe_df(df_dupl).write_excel(
            workbook=wb, worksheet="duplicate_reactions"
        )
        _excel_safe_df(df_imb).write_excel(workbook=wb, worksheet="curated_imbalances")
        df_ign.write_excel(workbook=wb, worksheet="ignored_imbalances")


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
    # cause 2): a faulty model can never freeze the run.
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
    Writes a PaRMMoSaHN curation DB format with the identified problems.
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

        return pl.DataFrame(df_rows)

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
            schema=["Pair", "Decision",
                    "ID 1", "Name 1", "Equation 1", "Lower bound 1", "Upper bound 1",
                    "ID 2", "Name 2", "Equation 2", "Lower bound 2", "Upper bound 2"]
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
) -> int:
    """Apply duplicate reaction removal decisions from the curation DB.

    Decision logic:
      0: Drop both reactions
      1: Keep ID 1, drop ID 2
      2: Keep ID 2, drop ID 1
      3: Keep both (do nothing)
      blank: Skip and continue

    Returns:
        Number of reactions removed.
    """
    changes_made = 0
    for row in df_dupl.iter_rows(named=True):
        decision = row.get("Decision")
        id1 = row.get("ID 1")
        id2 = row.get("ID 2")

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

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message = "need to pass in a list",
                category = UserWarning,
                )
            for rxn_id in to_drop:
                if rxn_id and rxn_id in model.reactions:
                    rxn_to_drop = cast(Reaction, model.reactions.get_by_id(rxn_id))
                    rxn_to_drop.remove_from_model(remove_orphans=True)
                    changes_made += 1

    return changes_made


def _apply_formula_fixes(
    model: "Model",
    df_cur: "pl.DataFrame",
    model_name: str,
) -> int:
    """Apply stoichiometric formula corrections from the curated_imbalances sheet.

    When *stats* is given, each row's hit count is accumulated under
    ``curated_imbalances`` keyed by ``reaction_id`` so the orchestrator can
    surface DB rows that fired on zero strains (i.e. stale corrections).

    Returns:
        Number of reactions updated.
    """
    changes_made = 0
    for row in df_cur.iter_rows(named=True):
        rxn_id = row.get("reaction_id")
        new_formula = row.get("new_formula")

        if rxn_id is None or new_formula is None or str(rxn_id) not in model.reactions:
            continue

        rxn = cast(Reaction, model.reactions.get_by_id(str(rxn_id)))
        old_formula = rxn.build_reaction_string()

        if str(new_formula).strip() != old_formula:
            try:
                rxn.build_reaction_from_string(str(new_formula).strip())
                changes_made += 1
            except Exception as e:
                logger.warning(
                    f"Could not apply new formula for {rxn_id} in {model_name}: {e}"
                )

    return changes_made


def apply_curations(
    models_dir: Path,
    curation_db_path: Path,
    output_dir: Path,
) -> Path | None:
    """Read the curation database and apply all fixes to a specific strain model.

    Orchestrates duplicate-reaction decisions and formula corrections.

    Skip-existing: if the output already exists it is returned untouched so a
    resumed run does not re-curate. The output is written atomically, so an
    existing file is always complete (never a truncated partial). To re-apply an
    updated curation database, clear the output directory first.
    """
    models = list(models_dir.rglob('*.xml'))
    for model_path in models:
        out_file = output_dir / model_path.name
        
        logger.debug(f"Curating model {out_file.stem}")
        
        try:
            model = cast(Model, cb.io.read_sbml_model(str(model_path)))
        except Exception as e:
            logger.error(f"Failed to load {model_path.name} for curation: {e}")
            return None
    
        try:
            df_dupl = pl.read_excel(
                curation_db_path, sheet_name="duplicate_reactions", engine="openpyxl"
            )
            df_cur = pl.read_excel(
                curation_db_path, sheet_name="curated_imbalances", engine="openpyxl"
            )
        except Exception as e:
            logger.error(f"Could not read curation database {curation_db_path}: {e}")
            return None
    
        changes_made = 0
        changes_made += _apply_duplicate_decisions(model, df_dupl)
        changes_made += _apply_formula_fixes(model, df_cur, model_path.name)
    
        # Atomic write: a crash mid-write must not leave a truncated SBML that the
        # skip-existing guard above would later mistake for a completed curation.
        tmp_out = out_file.with_suffix(out_file.suffix + ".tmp")
        cb.io.write_sbml_model(model, str(tmp_out))
        os.replace(tmp_out, out_file)
        logger.debug(f"Applied {changes_made} curations to {model_path.name}")
        
    logger.info(f"Curated {len(models)} models.")
    
    return output_dir
