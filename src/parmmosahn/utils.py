import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Generator, Optional

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
    If stdout_file or stderr_file are provided, outputs are written to them instead of captured.
    Returns standard output as string if stdout_file is None, else returns None.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    stdout_stream = None
    stderr_stream = None
    try:
        kwargs: dict[str, Any] = {"text": True, "check": True}

        # Determine stdout handling
        if stdout_file is not None:
            stdout_stream = stdout_file.open("a")
            kwargs["stdout"] = stdout_stream
        else:
            kwargs["stdout"] = subprocess.PIPE

        # Determine stderr handling
        if stderr_file is not None:
            stderr_stream = stderr_file.open("a")
            kwargs["stderr"] = stderr_stream
        else:
            kwargs["stderr"] = subprocess.PIPE

        result = subprocess.run(cmd, **kwargs)

        # Close files if opened
        if stdout_stream is not None:
            stdout_stream.close()
        if stderr_stream is not None:
            stderr_stream.close()

        # Return stdout if captured
        if stdout_file is None:
            return result.stdout
        return None

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        if stderr_file is None:
            logger.error(f"Standard Error output:\n{e.stderr}")
        else:
            logger.error(f"Check {stderr_file} for detailed error logs.")

        # Ensure we close file streams on error
        if stdout_stream is not None and not stdout_stream.closed:
            stdout_stream.close()
        if stderr_stream is not None and not stderr_stream.closed:
            stderr_stream.close()

        raise ToolExecutionError(f"Tool failed with exit code {e.returncode}") from e


def extract_proteome_from_gbff(gbff_path: Path, output_faa_path: Path) -> Path:
    """
    Extracts translated CDS features from a GenBank flat file into a FASTA file.

    Architectural notes:
    - Yields features generator-style to maintain low memory overhead.
    - Uses TemporaryDirectory for safe temporary file creation.
    """
    if output_faa_path.exists():
        logger.info(
            f"Proteome already exists at {output_faa_path}. Skipping extraction."
        )
        return output_faa_path

    def _generate_proteins() -> Generator[SeqRecord, None, None]:
        with open(gbff_path, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        translation = feature.qualifiers.get("translation", [None])[0]
                        if not translation:
                            continue

                        # Prefer protein_id, fallback to locus_tag, then feature ID
                        seq_id = feature.qualifiers.get(
                            "protein_id",
                            feature.qualifiers.get("locus_tag", ["unknown_id"]),
                        )[0]

                        yield SeqRecord(
                            Seq(translation),
                            id=seq_id,
                            description=f"Extracted from {record.id}",
                        )

    try:
        with tempfile.TemporaryDirectory(dir=output_faa_path.parent) as temp_dir:
            temp_path = Path(temp_dir) / "temp.faa"

            with open(temp_path, "w") as temp_out:
                count = SeqIO.write(_generate_proteins(), temp_out, "fasta")

            if count == 0:
                logger.warning(f"No CDS translations found in {gbff_path}!")

            # Move from temp dir to final target
            temp_path.replace(output_faa_path)
            logger.debug(
                f"Successfully extracted {count} proteins to {output_faa_path}"
            )

    except Exception as e:
        logger.error(f"Failed to extract proteome from {gbff_path}: {e}")
        raise

    return output_faa_path
