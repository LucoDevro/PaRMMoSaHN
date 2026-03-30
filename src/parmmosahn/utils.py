import logging
import subprocess
import tempfile
import threading
from pathlib import Path
from typing import Any, Generator, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class ToolExecutionError(Exception):
    pass


def run_command(
    cmd: list[str], stdout_file: Path | None = None, stderr_file: Path | None = None
) -> str | None:
    """
    Runs an external command robustly with REAL-TIME logging.
    Separates stdout and stderr completely to avoid stream pollution.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")

    out_f = stdout_file.open("a") if stdout_file else None
    err_f = stderr_file.open("a") if stderr_file else None

    captured_stdout: list[str] = []
    captured_stderr: list[str] = []

    def stream_reader(pipe, file_stream, log_func, capture_list):
        """Worker thread function to read streams without blocking."""
        for line in pipe:
            clean_line = line.strip()
            if clean_line:
                if file_stream is None:
                    log_func(f"[{cmd[0]}] {clean_line}")
                    capture_list.append(clean_line)
                else:
                    file_stream.write(line)
                    file_stream.flush()

    try:
        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True,
        ) as process:

            # Use threads to read both streams concurrently
            t_out = threading.Thread(
                target=stream_reader,
                args=(process.stdout, out_f, logger.info, captured_stdout),
            )
            t_err = threading.Thread(
                target=stream_reader,
                args=(process.stderr, err_f, logger.error, captured_stderr),
            )

            t_out.start()
            t_err.start()

            t_out.join()
            t_err.join()
            process.wait()

            if process.returncode != 0:
                # If it failed, bundle the stderr into the exception
                err_msg = (
                    "\n".join(captured_stderr)
                    if not stderr_file
                    else f"Check {stderr_file}."
                )
                raise subprocess.CalledProcessError(
                    process.returncode, cmd, stderr=err_msg
                )

        return "\n".join(captured_stdout) if stdout_file is None else None

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        if e.stderr:
            logger.error(f"Standard Error output:\n{e.stderr}")
        raise ToolExecutionError(f"Tool failed with exit code {e.returncode}") from e

    finally:
        if out_f:
            out_f.close()
        if err_f:
            err_f.close()


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
