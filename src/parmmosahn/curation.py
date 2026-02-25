import json
import logging
import re
from pathlib import Path
from typing import Any

import cobra as cb
import pandas as pd

from parmmosahn.utils import run_command

logger = logging.getLogger(__name__)

# ==========================================
# Phase 1: memote evaluation & DB Generation
# ==========================================


def run_memote_snapshot(model_path: Path, output_json: Path, output_html: Path) -> None:
    """Runs memote report snapshot on a single model."""
    logger.info(f"Running memote snapshot on {model_path.name}...")
    # Requires your modified memote version that supports --json
    cmd = [
        "memote",
        "report",
        "snapshot",
        "--json",
        str(output_json),
        "--filename",
        str(output_html),
        str(model_path),
    ]
    # We suppress stdout/stderr as memote is very noisy
    run_command(
        cmd,
        stdout_file=output_json.parent / "memote_run.log",
        stderr_file=output_json.parent / "memote_run.err",
    )


def extract_memote_problems(
    report_dir: Path, model_dir: Path, output_dir: Path
) -> None:
    """
    Parses memote JSON reports to find imbalances and duplicates.
    Replaces: make_curation_db.py
    """
    logger.info("Extracting problem balances and duplicates from memote reports...")
    output_dir.mkdir(parents=True, exist_ok=True)

    def read_test_results(test_key: str) -> dict:
        rxns: set[Any] = set()
        first_encountered = {}
        for report_file in report_dir.glob("*.json"):
            model_name = report_file.stem
            try:
                with report_file.open("r") as f:
                    cont = json.load(f)
                encountered_rxns = cont["tests"][test_key]["data"]

                if not encountered_rxns:
                    continue

                if isinstance(encountered_rxns[0], list):
                    encountered_rxns = [tuple(i) for i in encountered_rxns]

                new_rxns = set(encountered_rxns).difference(rxns)
                for r in new_rxns:
                    first_encountered[r] = model_name
                rxns = rxns.union(new_rxns)
            except (KeyError, IndexError, json.JSONDecodeError):
                logger.warning(f"Could not parse {test_key} in {report_file.name}")
        return first_encountered

    first_encounter_charge = read_test_results("test_reaction_charge_balance")
    first_encounter_mass = read_test_results("test_reaction_mass_balance")
    first_encounter_dupl = read_test_results("test_find_duplicate_reactions")

    def save_problem_balances(incidence_dict: dict, save_name: str):
        df = []
        for rxn_id, model_name in incidence_dict.items():
            try:
                network = cb.io.read_sbml_model(str(model_dir / f"{model_name}.xml"))
                rxn_obj = network.reactions.get_by_id(rxn_id)
                df.append(
                    {
                        "ID": rxn_obj.id,
                        "Name": rxn_obj.name,
                        "Equation": rxn_obj.build_reaction_string(),
                        "Imbalance": str(rxn_obj.check_mass_balance()).strip("{}"),
                        "First_encounter": model_name,
                    }
                )
            except Exception as e:
                logger.debug(f"Skipping {rxn_id} in {model_name}: {e}")

        if df:
            pd.DataFrame(df).to_csv(output_dir / save_name, sep="\t", index=False)

    def save_duplicates(incidence_dict: dict, save_name: str):
        df = []
        id_regex = re.compile(r"[0-9]{5}")
        for rxn_pair, model_name in incidence_dict.items():
            try:
                network = cb.io.read_sbml_model(str(model_dir / f"{model_name}.xml"))
                rxn_pair_obj = [network.reactions.get_by_id(r) for r in rxn_pair]

                # Extract numbers to find default reaction
                rxn_ids = [
                    int(id_regex.findall(r.id)[0]) if id_regex.findall(r.id) else 99999
                    for r in rxn_pair_obj
                ]
                min_idx = rxn_ids.index(min(rxn_ids))

                default_rxn = rxn_pair_obj.pop(min_idx)
                non_default_rxn = rxn_pair_obj[0]

                df.append(
                    {
                        "First_encounter": model_name,
                        "Pair": str(rxn_pair).strip("()"),
                        "ID 1": default_rxn.id,
                        "Name 1": default_rxn.name,
                        "Equation 1": default_rxn.build_reaction_string(),
                        "Lower bound 1": default_rxn.lower_bound,
                        "Upper bound 1": default_rxn.upper_bound,
                        "GPR 1": default_rxn.gene_reaction_rule,
                        "ID 2": non_default_rxn.id,
                        "Name 2": non_default_rxn.name,
                        "Equation 2": non_default_rxn.build_reaction_string(),
                        "Lower bound 2": non_default_rxn.lower_bound,
                        "Upper bound 2": non_default_rxn.upper_bound,
                        "GPR 2": non_default_rxn.gene_reaction_rule,
                    }
                )
            except Exception as e:
                logger.debug(f"Skipping pair {rxn_pair} in {model_name}: {e}")

        if df:
            pd.DataFrame(df).to_csv(output_dir / save_name, sep="\t", index=False)

    save_problem_balances(first_encounter_charge, "charge_imbalances.tsv")
    save_problem_balances(first_encounter_mass, "mass_imbalances.tsv")
    save_duplicates(first_encounter_dupl, "duplicate_reactions.tsv")


# ==========================================
# Phase 2: Applying the Curation Database
# ==========================================


def apply_curations(
    source_model_path: Path, curated_model_path: Path, curations_database: Path
) -> None:
    """
    Applies the Excel decision database to a model.
    Replaces: curate_models.py
    """
    logger.info(f"Applying curations to {source_model_path.name}...")

    try:
        imbalances = pd.read_excel(
            curations_database, sheet_name="corrected_balances"
        ).fillna("")
        duplicates = pd.read_excel(
            curations_database, sheet_name="duplicate_reactions"
        ).fillna("")
    except Exception as e:
        logger.error(f"Failed to read curation DB: {e}. Are the sheet names correct?")
        raise

    model = cb.io.read_sbml_model(str(source_model_path))

    # Correct Imbalances
    for rxn_name, form in zip(imbalances["reaction_id"], imbalances["new_formula"]):
        try:
            rxn = model.reactions.get_by_id(rxn_name)
            if form:
                logger.debug(f"Correcting imbalanced {rxn_name}")
                rxn.build_reaction_from_string(str(form))
            else:
                logger.debug(f"Removing imbalanced {rxn_name}")
                rxn.remove_from_model(remove_orphans=True)
        except KeyError:
            continue

    # Remove Duplicates
    for rxn_pair_names0, decision in zip(duplicates["Pair"], duplicates["Decision"]):
        try:
            # Safely parse the string tuple from the excel sheet e.g. "'rxn1', 'rxn2'"
            pair_str = str(rxn_pair_names0).replace("'", "").replace('"', "")
            rxn_pair_names = [r.strip() for r in pair_str.split(",")]

            if len(rxn_pair_names) != 2:
                continue

            rxn1 = model.reactions.get_by_id(rxn_pair_names[0])
            rxn2 = model.reactions.get_by_id(rxn_pair_names[1])

            logger.debug(
                f"Correcting duplicate pair ({rxn1.id}, {rxn2.id}) with decision {decision}"
            )

            if decision == 1:
                rxn2.remove_from_model(remove_orphans=True)
            elif decision == -1:
                rxn1.remove_from_model(remove_orphans=True)
            elif decision == 0:
                rxn1.remove_from_model(remove_orphans=True)
                rxn2.remove_from_model(remove_orphans=True)
            elif decision == 2:
                pass  # Keep both
            else:
                logger.warning(
                    f"Unknown decision {decision} for pair {rxn_pair_names0}"
                )

        except KeyError:
            continue

    # Save the curated model
    curated_model_path.parent.mkdir(parents=True, exist_ok=True)
    cb.io.write_sbml_model(model, str(curated_model_path))
    logger.info(f"Saved curated model to {curated_model_path}")
