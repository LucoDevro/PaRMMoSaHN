from pathlib import Path
from unittest.mock import MagicMock, patch

import cobra as cb
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


def test_apply_curations(dummy_sbml, dummy_curation_db, tmp_path):
    """
    Test that the apply_curations function correctly parses the Excel DB
    and applies stoichiometry fixes and duplicate removals to a real SBML model.
    """
    curated_model_path = tmp_path / "curated_model.xml"

    # Run the function using the mock fixtures!
    apply_curations(
        source_model_path=dummy_sbml,
        curated_model_path=curated_model_path,
        curations_database=dummy_curation_db,
    )

    # Assert the output file was created
    assert curated_model_path.exists()

    # Load the curated model to verify the internal logic worked
    curated_model = cb.io.read_sbml_model(str(curated_model_path))

    # 1. Check if the balance correction was applied to RXN_1
    rxn1 = curated_model.reactions.get_by_id("RXN_1")
    mbc = curated_model.metabolites.get_by_id("B_c")

    # Tell the type checker exactly what these objects are!
    assert isinstance(rxn1, cb.Reaction)
    assert isinstance(mbc, cb.Metabolite)

    # According to our mock DB, it should now have 2.0 of B_c
    assert rxn1.metabolites[mbc] == 2.0

    # 2. Check if the duplicate removal logic worked
    # RXN_DUPE_1 should still be there
    assert "RXN_DUPE_1" in curated_model.reactions

    # RXN_DUPE_2 should have been removed (because Decision == 1)
    assert "RXN_DUPE_2" not in curated_model.reactions
