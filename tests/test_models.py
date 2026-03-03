from pathlib import Path

import pytest

from parmmosahn.models import clean_fasta_file


def test_clean_fasta_file_with_bad_chars(tmp_path: Path):
    """Test that problematic characters are removed from FASTA files."""
    fasta_file = tmp_path / "test.faa"
    # Create a dummy fasta file with the problematic '’' character
    fasta_file.write_text(">seq1’\nATGC’C\n", encoding="utf-8")

    clean_fasta_file(fasta_file)

    # Check if it was cleaned
    cleaned_content = fasta_file.read_text(encoding="utf-8")
    assert "’" not in cleaned_content
    assert cleaned_content == ">seq1\nATGCC\n"


def test_clean_fasta_file_no_bad_chars(tmp_path: Path):
    """Test that valid FASTA files are left untouched."""
    fasta_file = tmp_path / "test.faa"
    original_content = ">seq1\nATGCC\n"
    fasta_file.write_text(original_content, encoding="utf-8")

    clean_fasta_file(fasta_file)

    assert fasta_file.read_text(encoding="utf-8") == original_content


def test_clean_fasta_file_decoding_error(tmp_path: Path):
    """Test that the function raises an error for non-text files."""
    binary_file = tmp_path / "test.faa"
    binary_file.write_bytes(b"\x80\x81\x82")  # Invalid UTF-8 bytes

    with pytest.raises(UnicodeDecodeError):
        clean_fasta_file(binary_file)
