import contextlib
import json
import logging
import os
import re
import subprocess
import time
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any, Dict, Tuple, cast

import cobra as cb
import pandas as pd
from cobra.core import Model, Reaction

# Import memote API directly
from memote.suite import api as memote_api
from memote.suite.reporting import ReportConfiguration

from parmmosahn.utils import ToolExecutionError, run_command

logger = logging.getLogger(__name__)

# ==========================================
# Phase 1: memote evaluation & DB Generation
# ==========================================


def run_memote_snapshot(model_path: Path, output_json: Path, output_html: Path) -> None:
    """
    Runs memote report snapshot on a single model via the native Python API.
    Saves both the JSON test results and the rendered HTML report.
    """
    logger.info(f"Running memote snapshot on {model_path.name}...")

    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    # 1. Validate and load model via memote's parser
    model_obj, sbml_ver, notifications = memote_api.validate_model(str(model_path))
    if model_obj is None:
        logger.error(
            f"Failed to load model {model_path.name} for memote due to SBML errors."
        )
        return

    # 2. Test model (suppress stdout to avoid pytest polluting the console)
    pytest_args = ["--tb", "no", "-q"]

    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            # Tell Pylance this returns a Tuple because results=True
            output = memote_api.test_model(
                model_obj, sbml_version=sbml_ver, results=True, pytest_args=pytest_args
            )
            # Safely unpack by casting to the expected Tuple signature
            code, results = cast(Tuple[int, Any], output)

    # 3. Generate and save the reports
    config = ReportConfiguration.load()

    # Save HTML report
    html_report = memote_api.snapshot_report(results, config=config, html=True)
    if html_report:
        with output_html.open("w", encoding="utf-8") as f:
            f.write(cast(str, html_report))

    # Save JSON report (test data)
    json_report = memote_api.snapshot_report(results, config=config, html=False)
    if json_report:
        with output_json.open("w", encoding="utf-8") as f:
            f.write(cast(str, json_report))

    logger.info(f"Memote snapshot completed for {model_path.name}.")


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

    # --- NEW: Cache loaded models so we don't load the same XML 50 times! ---
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

    def get_duplicates_df() -> pd.DataFrame:
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
                rxn_ids = [
                    int(match[0]) if (match := id_regex.findall(r.id)) else 99999
                    for r in rxn_pair_obj
                ]

                min_idx = rxn_ids.index(min(rxn_ids))
                r1 = cast(Reaction, rxn_pair_obj.pop(min_idx))
                r2 = cast(Reaction, rxn_pair_obj[0])

                df_rows.append(
                    {
                        "Pair": str(rxn_pair),
                        "Decision": "",  # Leave empty for human input
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
        return pd.DataFrame(df_rows)

    def get_imbalance_df() -> pd.DataFrame:
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
        return pd.DataFrame(df_rows)

    # Write DataFrames to Excel
    logger.info(f"Writing curation template to {output_excel.name}...")
    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        df_dupl = get_duplicates_df()
        if df_dupl.empty:
            df_dupl = pd.DataFrame(
                columns=["Pair", "Decision", "ID 1", "Name 1", "Equation 1"]
            )
        df_dupl.to_excel(writer, sheet_name="duplicate_reactions", index=False)

        df_imb = get_imbalance_df()
        if df_imb.empty:
            df_imb = pd.DataFrame(
                columns=[
                    "reaction_id",
                    "name",
                    "new_formula",
                    "balance_type",
                    "modification",
                    "Comment",
                ]
            )
        df_imb.to_excel(writer, sheet_name="curated_imbalances", index=False)

        pd.DataFrame(
            columns=[
                "reaction_id",
                "name",
                "Imbalance",
                "Reason to leave intact",
                "Comment",
            ]
        ).to_excel(writer, sheet_name="ignored_imbalances", index=False)


# ==========================================
# Phase 2: Applying the Curation Database
# ==========================================


def apply_curations(
    model_path: Path, curation_db_path: Path, output_dir: Path
) -> Path | None:
    """
    Reads the curation database and applies fixes to a specific strain model.
    Drops duplicate reactions based on 'Decision' and updates formulas based on 'new_formula'.
    """
    out_file = output_dir / model_path.name

    # Check if we've already curated this file to save time on pipeline restarts
    if out_file.exists():
        return out_file

    try:
        # Cast the loaded model so Pylance knows its type
        model = cast(Model, cb.io.read_sbml_model(str(model_path)))
    except Exception as e:
        logger.error(f"Failed to load {model_path.name} for curation: {e}")
        return None

    # Load Excel DB
    try:
        # Sheet 0 = Duplicates (Tab 1), Sheet 1 = Curations (Tab 2)
        df_dupl = pd.read_excel(curation_db_path, sheet_name=0)
        df_cur = pd.read_excel(curation_db_path, sheet_name=1)
    except Exception as e:
        logger.error(f"Could not read curation database {curation_db_path}: {e}")
        return None

    changes_made = 0

    # 1. Apply Duplicate Fixes
    for _, row in df_dupl.iterrows():
        decision = row.get("Decision")
        id1 = row.get("ID 1")
        id2 = row.get("ID 2")

        # Skip if human hasn't made a decision
        if pd.isna(decision) or not str(decision).strip():
            continue

        decision = str(decision).strip()
        to_drop = id2 if decision == "1" else id1 if decision == "2" else None

        if to_drop and to_drop in model.reactions:
            # Cast to Reaction before interacting with it
            rxn_to_drop = cast(Reaction, model.reactions.get_by_id(to_drop))
            rxn_to_drop.remove_from_model(remove_orphans=True)
            changes_made += 1

    # 2. Apply Formula Fixes
    for _, row in df_cur.iterrows():
        rxn_id = row.get("reaction_id")
        new_formula = row.get("new_formula")

        # Skip if no reaction ID, no formula, or reaction not in this specific strain model
        if (
            pd.isna(rxn_id)
            or pd.isna(new_formula)
            or str(rxn_id) not in model.reactions
        ):
            continue

        # Cast to Reaction so Pylance recognizes .build_reaction_string()
        rxn = cast(Reaction, model.reactions.get_by_id(str(rxn_id)))
        old_formula = rxn.build_reaction_string()

        # Only apply if the formula actually changed
        if str(new_formula).strip() != old_formula:
            try:
                # Cobrapy method to safely overwrite stoichiometry
                rxn.build_reaction_from_string(str(new_formula).strip())
                changes_made += 1
            except Exception as e:
                logger.warning(
                    f"Could not apply new formula for {rxn_id} in {model_path.name}: {e}"
                )

    # 3. Save the model
    cb.io.write_sbml_model(model, str(out_file))
    logger.debug(f"Applied {changes_made} curations to {model_path.name}")

    return out_file


def run_sbo_annotation(
    model_dir: str | Path,
    output_dir: str | Path,
    sbo_wrapper_script: str | Path,
    sbo_dbs_dir: str | Path,
) -> Path:
    """
    Runs the custom SBOannotator wrapper on all curated models.
    """
    model_dir = Path(model_dir).resolve()
    output_dir = Path(output_dir).resolve()
    sbo_wrapper = Path(sbo_wrapper_script).resolve()
    dbs_dir = Path(sbo_dbs_dir).resolve()

    if not sbo_wrapper.exists():
        raise FileNotFoundError(f"SBOannotator wrapper not found at {sbo_wrapper}")

    models_out = output_dir / "models"
    logs_out = output_dir / "logs"
    models_out.mkdir(parents=True, exist_ok=True)
    logs_out.mkdir(parents=True, exist_ok=True)
    dbs_dir.mkdir(parents=True, exist_ok=True)

    xml_models = list(model_dir.glob("*.xml"))
    logger.info(f"Adding SBO terms to {len(xml_models)} models...")

    for xml in xml_models:
        strain_name = xml.stem
        logger.info(f"Running SBOannotator for {strain_name}...")

        dest_xml = models_out / xml.name
        log_out = logs_out / f"{strain_name}.log"
        log_err = logs_out / f"{strain_name}.err"

        cmd = ["python", str(sbo_wrapper), str(xml), str(dest_xml), str(dbs_dir)]

        # Using our robust utils function to capture logs directly to files
        run_command(cmd, stdout_file=log_out, stderr_file=log_err)

    return models_out


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

    xml_models = list(model_dir.glob("*.xml"))
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

        # <-- NEW: Wait for PostgreSQL to actually finish booting up!
        logger.info("Waiting 15 seconds for databases to initialize...")
        time.sleep(240)

        # --- ADD THESE LINES TO FIX THE "FAST" ATTRIBUTE BUG ---
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
            for xml in xml_models:
                strain_name = xml.stem
                logger.info(f"Polishing {strain_name}...")

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
                    subprocess.run(
                        cmd,
                        cwd=str(compose_dir),
                        env=run_env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True,
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
    # Hide the raw reports in the strains folder to keep the root directory clean
    reports_dir = output_dir / "03-strains" / "reports_pre_curation"
    reports_dir.mkdir(parents=True, exist_ok=True)

    models = list(model_dir.rglob("*.xml"))
    logger.info(f"Running automated memote evaluation on {len(models)} draft models...")

    for model in models:
        out_json = reports_dir / f"{model.stem}.json"
        out_html = reports_dir / f"{model.stem}.html"
        run_memote_snapshot(model, out_json, out_html)

    out_excel = output_dir / "curation_template.xlsx"
    extract_memote_problems(reports_dir, model_dir, out_excel)
    return out_excel


def evaluate_final_models(model_dir: Path, output_dir: Path) -> None:
    """Hidden Evaluator: Runs memote on final models to prove FAIR compliance."""
    reports_dir = output_dir / "08-final_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    models = list(model_dir.rglob("*.xml"))
    if not models:
        logger.warning("No final XML models found to evaluate.")
        return

    logger.info(f"=== STEP 8: Final Memote Evaluation ===")
    logger.info(f"Running final memote evaluation on {len(models)} polished models...")
    for model in models:
        out_json = reports_dir / f"{model.stem}.json"
        out_html = reports_dir / f"{model.stem}.html"
        run_memote_snapshot(model, out_json, out_html)

    logger.info(f"Final FAIR-compliance reports saved to {reports_dir}")
