import contextlib
import json
import logging
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
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
from parmmosahn.utils import ToolExecutionError, run_command

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

    model_obj, sbml_ver, notifications = memote_api.validate_model(str(model_path))
    if model_obj is None:
        logger.error(
            f"Failed to load model {model_path.name} for memote due to SBML errors."
        )
        return

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
        with output_html.open("w", encoding="utf-8") as f:
            f.write(cast(str, html_report))

    json_report = memote_api.snapshot_report(results, config=config, html=False)
    if json_report:
        with output_json.open("w", encoding="utf-8") as f:
            f.write(cast(str, json_report))


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
        _run_memote_snapshot_inprocess(model_path, output_json, output_html)
        logger.debug(f"Memote snapshot completed for {model_path.name}.")
        return True

    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    with ProcessPoolExecutor(max_workers=1) as executor:
        future = executor.submit(
            _run_memote_snapshot_inprocess,
            model_path,
            output_json,
            output_html,
        )
        try:
            future.result(timeout=timeout)
            logger.debug(f"Memote snapshot completed for {model_path.name}.")
            return True
        except TimeoutError:
            logger.error(f"Memote timed out after {timeout}s for {model_path.name}.")
            return False
        except Exception as e:
            logger.error(f"Memote failed for {model_path.name}: {e}")
            return False


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
            except Exception:
                pass
        return first_encountered

    first_encounter_charge = read_test_results("test_reaction_charge_balance")
    first_encounter_mass = read_test_results("test_reaction_mass_balance")
    first_encounter_dupl = read_test_results("test_find_duplicate_reactions")

    # Trying to cache loaded models, too avoid loading the same XML 50 times.
    # TODO check if this actually improves performance or just adds complexity.
    _loaded_models = {}

    def get_cached_model(model_name: str):
        if model_name not in _loaded_models:
            # Check the new subfolder architecture
            model_path = model_dir / model_name / f"{model_name}.xml"
            # Fallback if it's flat
            if not model_path.exists():
                model_path = model_dir / f"{model_name}.xml"

            if model_path.exists():
                _loaded_models[model_name] = cast(
                    Model, cb.io.read_sbml_model(str(model_path))
                )
            else:
                _loaded_models[model_name] = None
        return _loaded_models[model_name]

    def get_duplicates_df() -> pl.DataFrame:
        df_rows = []
        id_regex = re.compile(r"\d{5}")
        for rxn_pair, model_name in first_encounter_dupl.items():
            network = get_cached_model(model_name)
            if not network:
                continue
            try:
                rxn_pair_obj = [
                    cast(Reaction, network.reactions.get_by_id(str(r)))
                    for r in rxn_pair
                ]

                # Supervisor fix: No 99999 fallback. Warn if pattern doesn't match!
                rxn_ids: list[float] = []
                for r in rxn_pair_obj:
                    match = id_regex.findall(r.id)
                    if match:
                        rxn_ids.append(float(match[0]))
                    else:
                        logger.warning(
                            f"Reaction ID '{r.id}' does not contain a 5-digit number for duplicate sorting."
                        )
                        rxn_ids.append(
                            float("inf")
                        )  # Ensure it sorts to the back safely without faking an ID

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
        df_dupl.write_excel(workbook=wb, worksheet="duplicate_reactions")
        df_imb.write_excel(workbook=wb, worksheet="curated_imbalances")
        df_ign.write_excel(workbook=wb, worksheet="ignored_imbalances")


# ==========================================
# Phase 2: Applying the Curation Database
# ==========================================


def _apply_duplicate_decisions(model: "Model", df_dupl: "pl.DataFrame") -> int:
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

        for rxn_id in to_drop:
            if rxn_id and rxn_id in model.reactions:
                rxn_to_drop = cast(Reaction, model.reactions.get_by_id(rxn_id))
                rxn_to_drop.remove_from_model(remove_orphans=True)
                changes_made += 1

    return changes_made


def _apply_formula_fixes(
    model: "Model", df_cur: "pl.DataFrame", model_name: str
) -> int:
    """Apply stoichiometric formula corrections from the curated_imbalances sheet.

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
    model_path: Path, curation_db_path: Path, output_dir: Path
) -> Path | None:
    """Read the curation database and apply all fixes to a specific strain model.

    Orchestrates duplicate-reaction decisions and formula corrections. The public
    signature is unchanged; internal logic is delegated to
    :func:`_apply_duplicate_decisions` and :func:`_apply_formula_fixes`.
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
        df_cur = pl.read_excel(
            curation_db_path, sheet_name="curated_imbalances", engine="openpyxl"
        )
    except Exception as e:
        logger.error(f"Could not read curation database {curation_db_path}: {e}")
        return None

    changes_made = 0
    changes_made += _apply_duplicate_decisions(model, df_dupl)
    changes_made += _apply_formula_fixes(model, df_cur, model_path.name)

    cb.io.write_sbml_model(model, str(out_file))
    logger.debug(f"Applied {changes_made} curations to {model_path.name}")

    return out_file


def run_model_polisher(
    model_dir: str | Path,
    output_dir: str | Path,
    engine: str = "docker",
) -> Path:
    """
    Executes the Palsson lab ModelPolisher via Docker Compose or alternative engines.
    Automatically provisions the local BiGG and Annotation databases.
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

        # Wait for PostgreSQL to actually finish booting up, otherwise it will reject connections from ModelPolisher and cause failures. 15 seconds was enough for newer versions, but 4 minutes seems necessary for older version.
        logger.info("Waiting 240 seconds for databases to initialize...")
        time.sleep(240)

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

                try:
                    result = subprocess.run(
                        cmd,
                        cwd=str(compose_dir),
                        env=run_env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True,
                    )
                    # Write per-model log even on success
                    log_file = logs_out / f"{strain_name}.log"
                    log_file.write_text(
                        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
                        encoding="utf-8",
                    )
                    success_count += 1

                except subprocess.CalledProcessError as e:
                    logger.error(f"ModelPolisher failed on {strain_name}.")
                    if e.stdout.strip():
                        logger.error(f"Docker STDOUT:\n{e.stdout.strip()}")
                    if e.stderr.strip():
                        logger.error(f"Docker STDERR:\n{e.stderr.strip()}")
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


def evaluate_final_models(model_dir: Path, output_dir: Path) -> None:
    """Hidden Evaluator: Runs memote on final models to prove FAIR compliance."""
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
    for i, model in enumerate(models, 1):
        logger.info(f"  [{i}/{n}] Evaluating {model.stem}...")
        out_json = reports_dir / f"{model.stem}.json"
        out_html = reports_dir / f"{model.stem}.html"
        if not run_memote_snapshot(model, out_json, out_html):
            failed.append(model.stem)

    if failed:
        logger.warning(f"Memote failed for {len(failed)}/{n} models: {failed}")
    logger.info(f"Final FAIR-compliance reports saved to {reports_dir}")
