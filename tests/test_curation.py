import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import cobra as cb
import pandas as pd
import pytest

from parmmosahn.curation import apply_curations, run_memote_snapshot


def test_run_memote_snapshot(tmp_path: Path):
    """
    Test that run_memote_snapshot successfully loads a model, suppresses output,
    generates reports via the memote API, and writes them to the correct paths.
    """
    # Setup paths
    model_path = tmp_path / "test_model.xml"
    output_json = tmp_path / "reports" / "test_model.json"
    output_html = tmp_path / "reports" / "test_model.html"

    # Create a dummy model file just so pathlib doesn't complain
    model_path.touch()

    # Create fake memote return data
    mock_model_obj = MagicMock()
    mock_sbml_ver = (3, 1, 2)
    mock_results = {"meta": "mocked test results"}

    # We use multiple patch decorators to intercept the heavy external dependencies
    with (
        patch("parmmosahn.curation.memote_api.validate_model") as mock_validate,
        patch("parmmosahn.curation.memote_api.test_model") as mock_test,
        patch("parmmosahn.curation.memote_api.snapshot_report") as mock_snapshot,
    ):

        # 1. Mock validate_model to return our fake tuple
        mock_validate.return_value = (mock_model_obj, mock_sbml_ver, None)

        # 2. Mock test_model to return an exit code (0) and our fake results dictionary
        mock_test.return_value = (0, mock_results)

        # 3. Mock snapshot_report to return HTML string if html=True, JSON string if html=False
        def side_effect_snapshot(results, config=None, html=True):
            return "<html>Mocked Report</html>" if html else '{"mock": "json"}'

        mock_snapshot.side_effect = side_effect_snapshot

        # Execute the function
        run_memote_snapshot(model_path, output_json, output_html)

        # Assertions
        mock_validate.assert_called_once_with(str(model_path))

        mock_test.assert_called_once_with(
            mock_model_obj,
            sbml_version=mock_sbml_ver,
            results=True,
            pytest_args=["--tb", "no", "-q"],
        )

        assert mock_snapshot.call_count == 2  # Once for HTML, once for JSON

        # Verify that the files were actually created and contain our mocked strings
        assert output_html.exists()
        assert output_json.exists()
        assert output_html.read_text(encoding="utf-8") == "<html>Mocked Report</html>"
        assert output_json.read_text(encoding="utf-8") == '{"mock": "json"}'


def test_run_memote_snapshot_invalid_model(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    """
    Test that the function gracefully aborts if memote fails to parse the SBML model.
    """
    model_path = tmp_path / "broken_model.xml"
    output_json = tmp_path / "broken_model.json"
    output_html = tmp_path / "broken_model.html"

    with patch("parmmosahn.curation.memote_api.validate_model") as mock_validate:
        # Simulate memote returning None for the model object (parsing failed)
        mock_validate.return_value = (None, None, {"errors": ["Invalid SBML"]})

        run_memote_snapshot(model_path, output_json, output_html)

        # Verify the error was logged and files were NOT created
        assert f"Failed to load model {model_path.name}" in caplog.text
        assert not output_json.exists()
        assert not output_html.exists()


@patch("parmmosahn.curation.cb.io.read_sbml_model")
def test_extract_memote_problems(mock_read_sbml, tmp_path):
    """
    Test that extract_memote_problems successfully parses JSON reports,
    loads the models to extract the formulas, and writes a correctly
    formatted 3-sheet Excel template.
    """
    from parmmosahn.curation import extract_memote_problems

    report_dir = tmp_path / "reports"
    model_dir = tmp_path / "models"
    report_dir.mkdir()
    model_dir.mkdir()

    fake_json = report_dir / "strain_A.json"
    fake_data = {
        "tests": {
            "test_reaction_charge_balance": {"data": ["RXN_C"]},
            "test_reaction_mass_balance": {"data": ["RXN_M"]},
            "test_find_duplicate_reactions": {"data": [["RXN_D1", "RXN_D2"]]},
        }
    }
    fake_json.write_text(json.dumps(fake_data), encoding="utf-8")

    # FIXED: Put the dummy XML inside the strain-specific subfolder!
    strain_dir = model_dir / "strain_A"
    strain_dir.mkdir()
    (strain_dir / "strain_A.xml").touch()

    mock_rxn = MagicMock()
    mock_rxn.id = "RXN_TEST"
    mock_rxn.name = "Fake Reaction"
    mock_rxn.build_reaction_string.return_value = "A -> B"
    mock_rxn.check_mass_balance.return_value = {"H": -1}
    mock_rxn.lower_bound = -1000
    mock_rxn.upper_bound = 1000
    mock_rxn.gene_reaction_rule = "(gene1)"

    mock_model = MagicMock()
    mock_model.reactions.get_by_id.return_value = mock_rxn
    mock_read_sbml.return_value = mock_model

    out_excel = tmp_path / "curation_template.xlsx"
    extract_memote_problems(report_dir, model_dir, out_excel)

    assert out_excel.exists()

    xls = pd.read_excel(out_excel, sheet_name=None)

    # Verify the sheets were created
    assert "duplicate_reactions" in xls
    assert "curated_imbalances" in xls
    assert "ignored_imbalances" in xls

    df_imbalances = xls["curated_imbalances"]
    assert len(df_imbalances) == 2  # 1 charge + 1 mass imbalance
    assert "balance_type" in df_imbalances.columns

    df_dupes = xls["duplicate_reactions"]
    assert len(df_dupes) == 1  # 1 duplicate pair


@patch("parmmosahn.curation.cb.io.read_sbml_model")
@patch("parmmosahn.curation.cb.io.write_sbml_model")
def test_apply_curations(mock_write, mock_read, tmp_path):
    """
    Test that apply_curations correctly parses the Excel DB
    and applies stoichiometry fixes and duplicate removals.
    """
    from parmmosahn.curation import apply_curations

    # 1. Setup mock paths
    model_path = tmp_path / "tiny_model.xml"
    model_path.touch()

    db_path = tmp_path / "curation_db.xlsx"
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    # 2. Create a fake Excel DB so pandas can read it
    with pd.ExcelWriter(db_path) as writer:
        pd.DataFrame(
            {"Decision": [1], "ID 1": ["RXN_D1"], "ID 2": ["RXN_D2"]}
        ).to_excel(writer, sheet_name="duplicate_reactions")

        pd.DataFrame(
            {"reaction_id": ["RXN_TEST"], "new_formula": ["A -> B + C"]}
        ).to_excel(writer, sheet_name="curated_imbalances")

    # 3. Setup mock cobrapy Model and Reactions
    mock_rxn_d1 = MagicMock()
    mock_rxn_d2 = MagicMock()
    mock_rxn_test = MagicMock()
    mock_rxn_test.build_reaction_string.return_value = "A -> B"

    mock_model = MagicMock()
    # Make the mock model say "Yes, I contain these reactions"
    mock_model.reactions.__contains__.side_effect = lambda x: x in [
        "RXN_D1",
        "RXN_D2",
        "RXN_TEST",
    ]

    def get_by_id_side_effect(rxn_id):
        if rxn_id == "RXN_D1":
            return mock_rxn_d1
        if rxn_id == "RXN_D2":
            return mock_rxn_d2
        if rxn_id == "RXN_TEST":
            return mock_rxn_test
        raise KeyError()

    mock_model.reactions.get_by_id.side_effect = get_by_id_side_effect
    mock_read.return_value = mock_model

    # 4. Execute the function
    result_path = apply_curations(
        model_path=model_path, curation_db_path=db_path, output_dir=out_dir
    )

    # 5. Assertions
    # Fixed: Verify the exact output path matches the input filename
    assert result_path == out_dir / "tiny_model.xml"

    # Verify duplicate was removed (Decision=1 means ID 2 is dropped)
    mock_rxn_d2.remove_from_model.assert_called_once_with(remove_orphans=True)

    # Verify formula was updated
    mock_rxn_test.build_reaction_from_string.assert_called_once_with("A -> B + C")

    # Verify the model was saved
    mock_write.assert_called_once_with(mock_model, str(result_path))


@patch("parmmosahn.curation.subprocess.run")
@patch("parmmosahn.curation.time.sleep")
@patch("parmmosahn.curation.files")
@patch("parmmosahn.curation.as_file")
def test_run_model_polisher(
    mock_as_file, mock_files, mock_sleep, mock_subprocess, tmp_path
):
    """
    Test that run_model_polisher correctly sets up the Docker Compose environment,
    waits for databases to initialize, and accurately loops over models.
    """
    from parmmosahn.curation import run_model_polisher

    model_dir = tmp_path / "draft_models"
    out_dir = tmp_path / "polished_models"
    model_dir.mkdir()

    # Create two dummy models
    (model_dir / "strain_A.xml").touch()
    (model_dir / "strain_B.xml").touch()

    # Mock the package data context manager for the docker-compose.yml file
    mock_compose_file = tmp_path / "docker-compose.yml"
    mock_compose_file.touch()
    mock_context = MagicMock()
    mock_context.__enter__.return_value = mock_compose_file
    mock_as_file.return_value = mock_context

    # Execute with podman engine to test engine branching
    results_path = run_model_polisher(model_dir, out_dir, engine="podman")

    # Verify outputs
    assert results_path == out_dir / "models"
    assert (out_dir / "models").exists()
    assert (out_dir / "logs").exists()

    # Verify sequence of events
    assert mock_sleep.called  # Should have waited for DB to boot

    # 1 'up' call, 2 'run' calls (for strain A and B), 1 'down' call = 4 total
    assert mock_subprocess.call_count == 4

    # Verify engine routing (should be using podman-compose)
    first_call_args = mock_subprocess.call_args_list[0][0][0]
    assert "podman-compose" in first_call_args
    assert "up" in first_call_args


def test_run_model_polisher_invalid_engine(tmp_path):
    """Test that the polisher rejects unsupported container engines."""
    from parmmosahn.curation import run_model_polisher

    with pytest.raises(ValueError, match="Unsupported container engine"):
        run_model_polisher(tmp_path, tmp_path, engine="virtualbox")


@patch("parmmosahn.curation.run_command")
def test_run_sbo_annotation(mock_run_command, tmp_path):
    """
    Test that the SBO annotation fallback correctly processes paths and executes.
    """
    from parmmosahn.curation import run_sbo_annotation

    model_dir = tmp_path / "in_models"
    out_dir = tmp_path / "sbo_out"
    model_dir.mkdir()

    (model_dir / "strain_X.xml").touch()

    sbo_wrapper = tmp_path / "SBOannotator.py"
    sbo_wrapper.touch()  # Must exist to pass the FileNotFoundError check
    sbo_dbs = tmp_path / "SBO_DBs"

    results_path = run_sbo_annotation(model_dir, out_dir, sbo_wrapper, sbo_dbs)

    assert results_path == out_dir / "models"
    assert (out_dir / "models").exists()
    assert (out_dir / "logs").exists()
    assert sbo_dbs.exists()

    # Verify that the wrapper function was called exactly once for our 1 model
    mock_run_command.assert_called_once()
