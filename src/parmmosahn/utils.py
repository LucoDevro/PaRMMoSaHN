import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Generator, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class ToolExecutionError(Exception):
    """Custom exception raised when an external bio-tool fails."""

    pass


def run_command(
    cmd: list[str], stdout_file: Path | None = None, stderr_file: Path | None = None
) -> str | None:
    """
    Runs an external command robustly.
    If stdout_file and stderr_file are provided, outputs are written to them instead of captured.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    try:
        # Check if BOTH files are provided for file-based output
        if stdout_file is not None and stderr_file is not None:
            with stdout_file.open("a") as out, stderr_file.open("a") as err:
                subprocess.run(cmd, stdout=out, stderr=err, text=True, check=True)
            return None

        # In all other cases (including if only one file is provided), capture to string
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        return result.stdout

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        if not (stdout_file and stderr_file):
            logger.error(f"Standard Error output:\n{e.stderr}")
        else:
            logger.error(f"Check {stderr_file} for detailed error logs.")
        raise ToolExecutionError(f"Tool failed with exit code {e.returncode}") from e


logger = logging.getLogger(__name__)


def extract_proteome_from_gbff(gbff_path: Path, output_faa_path: Path) -> Path:
    """
    Extracts translated CDS features from a GenBank flat file into a FASTA file.

    Architectural notes:
    - Yields features generator-style to maintain low memory overhead.
    - Atomic rename to prevent race conditions in HPC environments.
    """
    if output_faa_path.exists():
        logger.info(
            f"Proteome already exists at {output_faa_path}. Skipping extraction."
        )
        return output_faa_path

    # Define a generator to keep memory profile flat
    def _generate_proteins() -> Generator[SeqRecord, None, None]:
        with open(gbff_path, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        # Edge case: Pseudogenes or broken annotations might lack a translation
                        translation = feature.qualifiers.get("translation", [None])[0]
                        if not translation:
                            continue

                        # Edge case: Prefer locus_tag, fallback to protein_id, then feature ID
                        seq_id = feature.qualifiers.get(
                            "locus_tag",
                            feature.qualifiers.get("protein_id", ["unknown_id"]),
                        )[0]

                        yield SeqRecord(
                            Seq(translation),
                            id=seq_id,
                            description=f"Extracted from {record.id}",
                        )

    # Atomic write pattern for HPC safety
    fd, temp_path_str = tempfile.mkstemp(dir=output_faa_path.parent, suffix=".faa.tmp")
    temp_path = Path(temp_path_str)

    try:
        with os.fdopen(fd, "w") as temp_out:
            count = SeqIO.write(_generate_proteins(), temp_out, "fasta")

        if count == 0:
            logger.warning(f"No CDS translations found in {gbff_path}!")

        # Atomic rename to final destination
        temp_path.replace(output_faa_path)
        logger.debug(f"Successfully extracted {count} proteins to {output_faa_path}")

    except Exception as e:
        if temp_path.exists():
            temp_path.unlink()
        logger.error(f"Failed to extract proteome from {gbff_path}: {e}")
        raise

    return output_faa_path
