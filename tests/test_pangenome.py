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
    # (non-empty: build_pangenome rejects a zero-family FASTA).
    expected_extracted_fasta = output_dir / "softcore_protein_families.faa"
    expected_extracted_fasta.write_text(">fam1\nMAAA\n", encoding="utf-8")

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


@patch("parmmosahn.pangenome.run_command")
@patch("shutil.move")
def test_build_pangenome_soft_core_override(
    mock_move, mock_run_command, tmp_path: Path
):
    """soft_core_override should bypass the 2/n formula and use the given value."""
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    (genomes_dir / "test1.gbff").touch()
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    def mock_run_side_effect(cmd, **kwargs):
        if "info" in cmd:
            return "Genomes : 10\n"
        return ""

    mock_run_command.side_effect = mock_run_side_effect

    expected_fasta = output_dir / "softcore_protein_families.faa"
    expected_fasta.write_text(">fam1\nMAAA\n", encoding="utf-8")

    # Override with a specific threshold
    result_fasta = build_pangenome(
        genomes_dir=genomes_dir,
        output_dir=output_dir,
        label="override_test",
        threads=1,
        soft_core_override=0.75,
    )

    # fasta and write commands should still have run
    assert mock_run_command.call_count == 4

    # The fasta command should include the overridden threshold value
    fasta_call = [
        call for call in mock_run_command.call_args_list if "fasta" in call[0][0]
    ]
    assert len(fasta_call) == 1
    fasta_cmd = fasta_call[0][0][0]
    assert "0.75" in " ".join(str(a) for a in fasta_cmd)
