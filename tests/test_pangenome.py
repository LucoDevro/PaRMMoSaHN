from pathlib import Path
from unittest.mock import patch

import pytest

from parmmosahn.pangenome import build_pangenome, create_annotations_list


def test_create_annotations_list(tmp_path: Path):
    """Test that the annotations.list file is generated correctly with absolute paths."""
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()

    # Create dummy GenBank files
    gbff1 = genomes_dir / "genome1.gbff"
    gbff2 = genomes_dir / "genome2.gbff"
    gbff1.touch()
    gbff2.touch()

    output_list = tmp_path / "annotations.list"

    # Run the function
    create_annotations_list(genomes_dir, output_list)

    # Verify the file was created and contains the correct tab-separated data
    assert output_list.exists()
    content = output_list.read_text()
    assert "genome1.gbff" in content
    assert str(gbff1.resolve()) in content
    assert "genome2.gbff" in content


@patch("parmmosahn.pangenome.run_command")
@patch("shutil.move")
def test_build_pangenome_success(mock_move, mock_run_command, tmp_path: Path):
    """Test the full PPanGGOLiN workflow execution with mocked external tools."""
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    (
        genomes_dir / "test1.gbff"
    ).touch()  # Need at least one file so it doesn't fail early

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Create a smart mock function for run_command
    def mock_run_side_effect(cmd, **kwargs):
        # If the pipeline calls `ppanggolin info`, return fake terminal output
        if "info" in cmd:
            return "Some random text\nGenomes : 10\nOther data : 5\n"
        # For annotate, cluster, and fasta, just return empty string
        return ""

    mock_run_command.side_effect = mock_run_side_effect

    # Simulate PPanGGOLiN creating the expected soft-core fasta file
    expected_extracted_fasta = output_dir / "softcore_protein_families.faa"
    expected_extracted_fasta.touch()

    # Run the function
    result_fasta = build_pangenome(
        genomes_dir=genomes_dir, output_dir=output_dir, label="test_pan", threads=2
    )

    # Assertions
    assert mock_run_command.call_count == 4  # annotate, cluster, info, fasta

    # Verify that shutil.move was called to rename the file to "test_pan.faa"
    mock_move.assert_called_once_with(
        expected_extracted_fasta, output_dir / "test_pan.faa"
    )

    # Verify the final returned path is correct
    assert result_fasta == output_dir / "test_pan.faa"
