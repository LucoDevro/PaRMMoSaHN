from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from parmmosahn.cli import cli


@pytest.fixture
def runner():
    """Fixture to provide a Click CliRunner for all tests."""
    return CliRunner()


@patch("parmmosahn.cli.build_pangenome")
def test_cli_pangenome(mock_build_pangenome, runner, tmp_path: Path):
    """Test the 'pangenome' CLI command."""
    # The CLI requires genomes-dir to exist, so we create a dummy one
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()

    # Simulate the user typing: parmmosahn pangenome --genomes-dir ... --label test_label
    result = runner.invoke(
        cli,
        [
            "pangenome",
            "--genomes-dir",
            str(genomes_dir),
            "--label",
            "test_label",
            "--threads",
            "4",
        ],
    )

    # Assert the command succeeded without crashing
    assert result.exit_code == 0
    assert "Building pangenome 'test_label'..." in result.output

    # Assert the underlying Python function was called with the parsed arguments
    mock_build_pangenome.assert_called_once_with(
        str(genomes_dir), "01-pangenome", "test_label", threads=4
    )


def test_cli_pangenome_missing_args(runner):
    """Test that the CLI correctly rejects commands missing required arguments."""
    # Simulate missing the required --label and --genomes-dir flags
    result = runner.invoke(cli, ["pangenome"])

    assert result.exit_code != 0
    assert "Missing option '--genomes-dir'" in result.output


@patch("parmmosahn.pipeline.run_all")
def test_cli_run_all(mock_run_all, runner, tmp_path: Path):
    """Test the 'run' command (the main pipeline orchestrator)."""
    # Create required dummy directories and files
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    medium_file = tmp_path / "medium.csv"
    medium_file.touch()
    db_file = tmp_path / "curation_db.xlsx"
    db_file.touch()

    output_dir = tmp_path / "results"

    # Simulate the user executing the full pipeline
    result = runner.invoke(
        cli,
        [
            "run",
            "--genomes-dir",
            str(genomes_dir),
            "--output-dir",
            str(output_dir),
            "--label",
            "my_strains",
            "--medium",
            str(medium_file),
            "--curation-db",
            str(db_file),
            "--threads",
            "8",
            "--skip-polisher",
        ],
    )

    assert result.exit_code == 0, f"CLI command failed: {result.output}"

    # FIXED: The CLI passes arguments as explicitly named keyword arguments (kwargs)
    mock_run_all.assert_called_once_with(
        genomes_dir=Path(genomes_dir),
        output_dir=Path(output_dir),
        label="my_strains",
        medium_file=Path(medium_file),
        curation_db=Path(db_file),
        threads=8,
        skip_polisher=True,
        engine="docker",  # The default value in the click option
    )
