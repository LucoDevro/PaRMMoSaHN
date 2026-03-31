from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from parmmosahn.results import gather_and_convert_models


def test_gather_finds_nested_models(tmp_path: Path):
    """gather_and_convert_models must find .xml files in strain subdirectories (rglob fix)."""
    model_dir = tmp_path / "models"
    out_dir = tmp_path / "final"

    # Create nested layout: models/strain_A/strain_A.xml
    strain_a = model_dir / "strain_A"
    strain_a.mkdir(parents=True)
    (strain_a / "strain_A.xml").touch()

    strain_b = model_dir / "strain_B"
    strain_b.mkdir()
    (strain_b / "strain_B.xml").touch()

    mock_model = MagicMock()

    with patch("parmmosahn.results.cb.io.read_sbml_model", return_value=mock_model):
        with patch("parmmosahn.results.cb.io.save_json_model"):
            with patch("parmmosahn.results.cb.io.save_matlab_model"):
                result = gather_and_convert_models(model_dir, out_dir)

    assert result == out_dir
    # Both XML files should have been copied to output_dir
    assert (out_dir / "strain_A.xml").exists()
    assert (out_dir / "strain_B.xml").exists()


def test_gather_flat_layout_still_works(tmp_path: Path):
    """gather_and_convert_models should also handle a flat directory layout."""
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    out_dir = tmp_path / "final"
    (model_dir / "polished.xml").touch()

    mock_model = MagicMock()
    with patch("parmmosahn.results.cb.io.read_sbml_model", return_value=mock_model):
        with patch("parmmosahn.results.cb.io.save_json_model"):
            with patch("parmmosahn.results.cb.io.save_matlab_model"):
                result = gather_and_convert_models(model_dir, out_dir)

    assert (out_dir / "polished.xml").exists()


def test_gather_empty_directory_logs_warning(tmp_path: Path, caplog):
    """gather_and_convert_models should warn and return early when no models are found."""
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    out_dir = tmp_path / "final"

    result = gather_and_convert_models(model_dir, out_dir)

    assert result == out_dir
    assert "No .xml models found" in caplog.text


def test_gather_handles_conversion_failure(tmp_path: Path, caplog):
    """A model that fails cobrapy loading should be logged but not crash the run."""
    model_dir = tmp_path / "models"
    out_dir = tmp_path / "final"
    strain_dir = model_dir / "bad"
    strain_dir.mkdir(parents=True)
    (strain_dir / "bad.xml").touch()

    with patch(
        "parmmosahn.results.cb.io.read_sbml_model",
        side_effect=Exception("parse error"),
    ):
        result = gather_and_convert_models(model_dir, out_dir)

    assert result == out_dir
    assert "Failed to convert" in caplog.text
