import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import cobra as cb
import pandas as pd
import pytest

from parmmosahn.curation import (
    _apply_duplicate_decisions,
    _apply_formula_fixes,
    apply_curations,
    run_memote_snapshot,
)

# ---------------------------------------------------------------------------
# run_memote_snapshot
# ---------------------------------------------------------------------------


def test_run_memote_snapshot(tmp_path: Path):
    """run_memote_snapshot should load a model, call memote API, and write reports."""
    model_path = tmp_path / "test_model.xml"
    output_json = tmp_path / "reports" / "test_model.json"
    output_html = tmp_path / "reports" / "test_model.html"
    model_path.touch()

    mock_model_obj = MagicMock()
    mock_sbml_ver = (3, 1, 2)
    mock_results = {"meta": "mocked test results"}

    with (
        patch("parmmosahn.curation.memote_api.validate_model") as mock_validate,
        patch("parmmosahn.curation.memote_api.test_model") as mock_test,
        patch("parmmosahn.curation.memote_api.snapshot_report") as mock_snapshot,
    ):
        mock_validate.return_value = (mock_model_obj, mock_sbml_ver, None)
        mock_test.return_value = (0, mock_results)

        def side_effect_snapshot(results, config=None, html=True):
            return "<html>Mocked Report</html>" if html else '{"mock": "json"}'

        mock_snapshot.side_effect = side_effect_snapshot

        run_memote_snapshot(model_path, output_json, output_html)

        mock_validate.assert_called_once_with(str(model_path))
        mock_test.assert_called_once_with(
            mock_model_obj,
            sbml_version=mock_sbml_ver,
            results=True,
            pytest_args=["--tb", "no", "-q"],
        )
        assert mock_snapshot.call_count == 2

        assert output_html.exists()
        assert output_json.exists()
        assert output_html.read_text(encoding="utf-8") == "<html>Mocked Report</html>"
        assert output_json.read_text(encoding="utf-8") == '{"mock": "json"}'


def test_run_memote_snapshot_invalid_model(tmp_path: Path, caplog):
    """run_memote_snapshot should abort gracefully when memote cannot parse the model."""
    model_path = tmp_path / "broken_model.xml"
    output_json = tmp_path / "broken_model.json"
    output_html = tmp_path / "broken_model.html"

    with patch("parmmosahn.curation.memote_api.validate_model") as mock_validate:
        mock_validate.return_value = (None, None, {"errors": ["Invalid SBML"]})
        run_memote_snapshot(model_path, output_json, output_html)

        assert f"Failed to load model {model_path.name}" in caplog.text
        assert not output_json.exists()
        assert not output_html.exists()


# ---------------------------------------------------------------------------
# extract_memote_problems
# ---------------------------------------------------------------------------


@patch("parmmosahn.curation.cb.io.read_sbml_model")
def test_extract_memote_problems(mock_read_sbml, tmp_path: Path):
    """extract_memote_problems should parse JSON reports and write a 3-sheet Excel.

    The polars write_excel call requires xlsxwriter (optional dep). We mock it so
    the test focuses on the data-extraction logic, not on the Excel writer.
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

    # Capture the DataFrames that would be written to Excel without needing xlsxwriter
    written_sheets: dict = {}

    class _FakeWorkbook:
        """Minimal workbook stub that records written DataFrames by sheet name."""

        def close(self):
            pass

    def fake_write_excel(df, workbook, worksheet, **kwargs):
        if isinstance(workbook, Path):
            written_sheets[worksheet] = df
            return _FakeWorkbook()
        elif isinstance(workbook, _FakeWorkbook):
            written_sheets[worksheet] = df
            return workbook

    import polars as pl

    with patch.object(pl.DataFrame, "write_excel", fake_write_excel):
        out_excel = tmp_path / "curation_template.xlsx"
        extract_memote_problems(report_dir, model_dir, out_excel)

    assert "duplicate_reactions" in written_sheets
    assert "curated_imbalances" in written_sheets
    assert "ignored_imbalances" in written_sheets
    assert len(written_sheets["curated_imbalances"]) == 2  # 1 charge + 1 mass
    assert len(written_sheets["duplicate_reactions"]) == 1


# ---------------------------------------------------------------------------
# _apply_duplicate_decisions (unit test of private helper)
# ---------------------------------------------------------------------------


def _make_mock_model(*reaction_ids):
    """Build a minimal mock cobra model that contains the given reaction IDs."""
    mock_model = MagicMock()
    reactions = {}
    for rxn_id in reaction_ids:
        rxn = MagicMock()
        rxn.id = rxn_id
        reactions[rxn_id] = rxn
    mock_model.reactions.__contains__.side_effect = lambda x: x in reactions
    mock_model.reactions.get_by_id.side_effect = lambda x: reactions[x]
    return mock_model, reactions


def test_apply_duplicate_decisions_decision_1():
    """Decision '1': keep ID 1, drop ID 2."""
    import polars as pl

    model, rxns = _make_mock_model("RXN_A", "RXN_B")
    df = pl.DataFrame({"Decision": ["1"], "ID 1": ["RXN_A"], "ID 2": ["RXN_B"]})

    changes = _apply_duplicate_decisions(model, df)

    assert changes == 1
    rxns["RXN_A"].remove_from_model.assert_not_called()
    rxns["RXN_B"].remove_from_model.assert_called_once_with(remove_orphans=True)


def test_apply_duplicate_decisions_decision_2():
    """Decision '2': keep ID 2, drop ID 1."""
    import polars as pl

    model, rxns = _make_mock_model("RXN_A", "RXN_B")
    df = pl.DataFrame({"Decision": ["2"], "ID 1": ["RXN_A"], "ID 2": ["RXN_B"]})

    changes = _apply_duplicate_decisions(model, df)

    assert changes == 1
    rxns["RXN_A"].remove_from_model.assert_called_once_with(remove_orphans=True)
    rxns["RXN_B"].remove_from_model.assert_not_called()


def test_apply_duplicate_decisions_decision_0():
    """Decision '0': drop both reactions."""
    import polars as pl

    model, rxns = _make_mock_model("RXN_A", "RXN_B")
    df = pl.DataFrame({"Decision": ["0"], "ID 1": ["RXN_A"], "ID 2": ["RXN_B"]})

    changes = _apply_duplicate_decisions(model, df)

    assert changes == 2
    rxns["RXN_A"].remove_from_model.assert_called_once_with(remove_orphans=True)
    rxns["RXN_B"].remove_from_model.assert_called_once_with(remove_orphans=True)


def test_apply_duplicate_decisions_decision_3():
    """Decision '3': keep both, no changes."""
    import polars as pl

    model, rxns = _make_mock_model("RXN_A", "RXN_B")
    df = pl.DataFrame({"Decision": ["3"], "ID 1": ["RXN_A"], "ID 2": ["RXN_B"]})

    changes = _apply_duplicate_decisions(model, df)

    assert changes == 0
    rxns["RXN_A"].remove_from_model.assert_not_called()
    rxns["RXN_B"].remove_from_model.assert_not_called()


def test_apply_duplicate_decisions_blank_skipped():
    """Blank decision should be silently skipped."""
    import polars as pl

    model, rxns = _make_mock_model("RXN_A", "RXN_B")
    df = pl.DataFrame({"Decision": [None], "ID 1": ["RXN_A"], "ID 2": ["RXN_B"]})

    changes = _apply_duplicate_decisions(model, df)

    assert changes == 0


# ---------------------------------------------------------------------------
# _apply_formula_fixes (unit test of private helper)
# ---------------------------------------------------------------------------


def test_apply_formula_fixes_updates_formula():
    """_apply_formula_fixes should call build_reaction_from_string when formula differs."""
    import polars as pl

    model, rxns = _make_mock_model("RXN_TEST")
    rxns["RXN_TEST"].build_reaction_string.return_value = "A -> B"
    df = pl.DataFrame({"reaction_id": ["RXN_TEST"], "new_formula": ["A -> B + C"]})

    changes = _apply_formula_fixes(model, df, "test_model.xml")

    assert changes == 1
    rxns["RXN_TEST"].build_reaction_from_string.assert_called_once_with("A -> B + C")


def test_apply_formula_fixes_skips_unchanged():
    """_apply_formula_fixes should skip a reaction whose formula is already correct."""
    import polars as pl

    model, rxns = _make_mock_model("RXN_TEST")
    rxns["RXN_TEST"].build_reaction_string.return_value = "A -> B"
    df = pl.DataFrame({"reaction_id": ["RXN_TEST"], "new_formula": ["A -> B"]})

    changes = _apply_formula_fixes(model, df, "test_model.xml")

    assert changes == 0
    rxns["RXN_TEST"].build_reaction_from_string.assert_not_called()


def test_apply_formula_fixes_skips_missing_reaction():
    """_apply_formula_fixes should silently skip reactions not in the model."""
    import polars as pl

    model, _ = _make_mock_model()  # model has no reactions
    df = pl.DataFrame({"reaction_id": ["GHOST_RXN"], "new_formula": ["A -> B"]})

    changes = _apply_formula_fixes(model, df, "test_model.xml")
    assert changes == 0


# ---------------------------------------------------------------------------
# apply_curations (integration: orchestrates both helpers)
# ---------------------------------------------------------------------------


@patch("parmmosahn.curation.cb.io.read_sbml_model")
@patch("parmmosahn.curation.cb.io.write_sbml_model")
def test_apply_curations(mock_write, mock_read, tmp_path: Path):
    """apply_curations should apply both duplicate decisions and formula fixes."""
    import polars as pl

    model_path = tmp_path / "tiny_model.xml"
    model_path.touch()
    db_path = tmp_path / "curation_db.xlsx"
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    mock_rxn_d1 = MagicMock()
    mock_rxn_d2 = MagicMock()
    mock_rxn_test = MagicMock()
    mock_rxn_test.build_reaction_string.return_value = "A -> B"

    mock_model = MagicMock()
    mock_model.reactions.__contains__.side_effect = lambda x: x in [
        "RXN_D1",
        "RXN_D2",
        "RXN_TEST",
    ]

    def get_by_id(rxn_id):
        return {
            "RXN_D1": mock_rxn_d1,
            "RXN_D2": mock_rxn_d2,
            "RXN_TEST": mock_rxn_test,
        }[rxn_id]

    mock_model.reactions.get_by_id.side_effect = get_by_id
    mock_read.return_value = mock_model

    df_dupl = pl.DataFrame({"Decision": ["1"], "ID 1": ["RXN_D1"], "ID 2": ["RXN_D2"]})
    df_cur = pl.DataFrame({"reaction_id": ["RXN_TEST"], "new_formula": ["A -> B + C"]})

    # Mock pl.read_excel — avoids the fastexcel / xlsxwriter optional deps
    with patch("parmmosahn.curation.pl.read_excel", side_effect=[df_dupl, df_cur]):
        result_path = apply_curations(
            model_path=model_path, curation_db_path=db_path, output_dir=out_dir
        )

    assert result_path == out_dir / "tiny_model.xml"
    mock_rxn_d2.remove_from_model.assert_called_once_with(remove_orphans=True)
    mock_rxn_test.build_reaction_from_string.assert_called_once_with("A -> B + C")
    mock_write.assert_called_once_with(mock_model, str(result_path))


@patch("parmmosahn.curation.cb.io.read_sbml_model")
@patch("parmmosahn.curation.cb.io.write_sbml_model")
def test_apply_curations_empty_curation_db(mock_write, mock_read, tmp_path: Path):
    """apply_curations should pass model through unchanged when curation DB has empty sheets."""
    import polars as pl

    model_path = tmp_path / "model.xml"
    model_path.touch()
    db_path = tmp_path / "empty_db.xlsx"
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    mock_model = MagicMock()
    mock_read.return_value = mock_model

    # Empty DataFrames with correct column schemas
    df_dupl = pl.DataFrame({"Decision": [], "ID 1": [], "ID 2": []}).cast(
        {"Decision": pl.Utf8, "ID 1": pl.Utf8, "ID 2": pl.Utf8}
    )
    df_cur = pl.DataFrame({"reaction_id": [], "new_formula": []}).cast(
        {"reaction_id": pl.Utf8, "new_formula": pl.Utf8}
    )

    with patch("parmmosahn.curation.pl.read_excel", side_effect=[df_dupl, df_cur]):
        result_path = apply_curations(
            model_path=model_path, curation_db_path=db_path, output_dir=out_dir
        )

    assert result_path == out_dir / "model.xml"
    mock_write.assert_called_once_with(mock_model, str(result_path))
    # No reactions should have been modified
    mock_model.reactions.get_by_id.assert_not_called()


def test_apply_curations_skips_existing_output(tmp_path: Path):
    """apply_curations should return early if the output file already exists."""
    model_path = tmp_path / "model.xml"
    model_path.touch()
    db_path = tmp_path / "db.xlsx"
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    # Pre-create the output file
    (out_dir / "model.xml").touch()

    with patch("parmmosahn.curation.cb.io.read_sbml_model") as mock_read:
        result = apply_curations(model_path, db_path, out_dir)

    assert result == out_dir / "model.xml"
    mock_read.assert_not_called()


# ---------------------------------------------------------------------------
# run_model_polisher
# ---------------------------------------------------------------------------


@patch("parmmosahn.curation.subprocess.run")
@patch("parmmosahn.curation.time.sleep")
@patch("parmmosahn.curation.files")
@patch("parmmosahn.curation.as_file")
def test_run_model_polisher(
    mock_as_file, mock_files, mock_sleep, mock_subprocess, tmp_path: Path
):
    """run_model_polisher should set up the container environment and loop over models."""
    from parmmosahn.curation import run_model_polisher

    model_dir = tmp_path / "draft_models"
    out_dir = tmp_path / "polished_models"
    model_dir.mkdir()
    (model_dir / "strain_A.xml").touch()
    (model_dir / "strain_B.xml").touch()

    mock_compose_file = tmp_path / "docker-compose.yml"
    mock_compose_file.touch()
    mock_context = MagicMock()
    mock_context.__enter__.return_value = mock_compose_file
    mock_as_file.return_value = mock_context

    results_path = run_model_polisher(model_dir, out_dir, engine="podman")

    assert results_path == out_dir / "models"
    assert (out_dir / "models").exists()
    assert (out_dir / "logs").exists()
    assert mock_sleep.called
    # 1 'up' + 2 'run' (strain A + B) + 1 'down' = 4 calls
    assert mock_subprocess.call_count == 4
    first_cmd = mock_subprocess.call_args_list[0][0][0]
    assert "podman-compose" in first_cmd
    assert "up" in first_cmd


def test_run_model_polisher_invalid_engine(tmp_path: Path):
    """run_model_polisher should reject unsupported container engines."""
    from parmmosahn.curation import run_model_polisher

    with pytest.raises(ValueError, match="Unsupported container engine"):
        run_model_polisher(tmp_path, tmp_path, engine="virtualbox")
