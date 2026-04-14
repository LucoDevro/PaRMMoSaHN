import contextlib
import logging
import shutil
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

    # 1. Start the ExitStack
    with contextlib.ExitStack() as stack:
        out_f = stack.enter_context(stdout_file.open("a")) if stdout_file else None
        err_f = stack.enter_context(stderr_file.open("a")) if stderr_file else None

        # 2. INDENT EVERYTHING INSIDE THE EXITSTACK
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
                    args=(process.stderr, err_f, logger.debug, captured_stderr),
                )

                t_out.start()
                t_err.start()

                # Wait for process first, THEN wait for threads to finish reading
                process.wait()
                t_out.join()
                t_err.join()

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

            error_details = [f"Tool failed with exit code {e.returncode}"]
            if e.stderr:
                error_details.append(f"stderr: {e.stderr.strip()}")
            elif stderr_file is not None:
                error_details.append(f"see stderr log at {stderr_file}")

            raise ToolExecutionError("; ".join(error_details)) from e


def validate_input_genomes(genomes_dir: Path) -> None:
    """Pre-flight validation of input genome files.

    Warns about potential issues that would cause pipeline failures later.
    Raises FileNotFoundError if no genome files found.
    """
    gbff_files = list(genomes_dir.glob("*.gbff")) + list(genomes_dir.glob("*.gbk"))
    faa_files = list(genomes_dir.glob("*.faa"))

    all_files = gbff_files + faa_files
    if not all_files:
        raise FileNotFoundError(
            f"No genome files (.gbff, .gbk, .faa) found in {genomes_dir}"
        )

    n_genomes = len(set(f.stem for f in all_files))
    if n_genomes < 3:
        logger.warning(
            f"Only {n_genomes} genome(s) found. Pangenome analysis typically "
            f"requires at least 3 genomes to produce meaningful results."
        )

    # Spot-check GenBank files for CDS features with translations
    for gbff in gbff_files[:3]:  # check up to 3 files
        cds_count = 0
        translated_count = 0
        try:
            for record in SeqIO.parse(str(gbff), "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        cds_count += 1
                        if feature.qualifiers.get("translation"):
                            translated_count += 1
                # Only scan first record per file for speed
                break
        except Exception as e:
            logger.warning(f"Could not parse {gbff.name}: {e}")
            continue

        if cds_count == 0:
            logger.warning(
                f"{gbff.name}: no CDS features found. This file may not be a "
                f"valid annotated genome (GenBank flat file format expected)."
            )
        elif translated_count == 0:
            logger.warning(
                f"{gbff.name}: CDS features found but none have translations. "
                f"Proteome extraction will produce empty output."
            )
        elif translated_count < 100:
            logger.warning(
                f"{gbff.name}: only {translated_count} translated CDS features "
                f"in first record. This may indicate a partial or low-quality genome."
            )


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


def extract_gene_id_mapping(gbff_path: Path, output_tsv: Path) -> Path:
    """Build a protein_id → NCBI RefSeq accession mapping from a GenBank file.

    Parses CDS features and extracts ``db_xref`` entries that start with
    ``RefSeq:``.  Writes a two-column TSV (protein_id \\t refseq_accession).
    """
    mapping: dict[str, str] = {}
    with open(gbff_path, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                protein_id = feature.qualifiers.get("protein_id", [None])[0]
                if not protein_id:
                    continue
                for xref in feature.qualifiers.get("db_xref", []):
                    if xref.startswith("RefSeq:"):
                        mapping[protein_id] = xref.split(":", 1)[1]
                        break

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.write_text(
        "\n".join(f"{k}\t{v}" for k, v in mapping.items()) + "\n",
        encoding="utf-8",
    )
    logger.debug(
        f"Gene ID mapping: {len(mapping)} RefSeq accessions from {gbff_path.name}"
    )
    return output_tsv


def load_config(config_path: Path) -> dict:
    """Load a YAML configuration file.

    Requires PyYAML (install with ``pip install 'parmmosahn[config]'``).

    Args:
        config_path: Path to the YAML config file.

    Returns:
        Dictionary of configuration values. Empty dict if the file is empty.

    Raises:
        ImportError: If PyYAML is not installed.
    """
    try:
        import yaml  # type: ignore[import]
    except ImportError as exc:
        raise ImportError(
            "PyYAML is required for config file support. "
            "Install it with: pip install 'parmmosahn[config]'  or  pip install pyyaml"
        ) from exc

    with config_path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    return data if isinstance(data, dict) else {}


def _extract_version(output: str) -> str:
    """Return the first line that looks like a version string, else the first non-empty line."""
    import re

    lines = [l.strip() for l in output.splitlines() if l.strip()]
    # Prefer a line that contains a semver-like pattern (e.g. '2.2.6', 'v1.4.2')
    version_pattern = re.compile(r"\bv?\d+\.\d+")
    for line in lines:
        if version_pattern.search(line):
            return line
    return lines[0] if lines else "unknown"


def check_environment() -> dict[str, dict[str, str | bool | None]]:
    """Check for required external tools and container engines.

    Probes each tool with its ``--version`` flag (30 s timeout to accommodate
    Python-based tools like PPanGGOLiN that have slow startup import times).

    Returns:
        Nested dict of ``{tool_name: {found: bool, version: str | None, path: str | None}}``.
    """
    tools: dict[str, list[str]] = {
        "ppanggolin": ["ppanggolin", "--version"],
        "gapseq": ["gapseq", "-v"],
        "diamond": ["diamond", "--version"],
        "docker": ["docker", "--version"],
        "podman": ["podman", "--version"],
        "apptainer": ["apptainer", "--version"],
    }

    results: dict[str, dict[str, str | bool | None]] = {}

    for name, version_cmd in tools.items():
        path = shutil.which(name)
        if path is None:
            results[name] = {"found": False, "version": None, "path": None}
            continue
        try:
            proc = subprocess.run(
                version_cmd,
                capture_output=True,
                text=True,
                timeout=30,
            )
            combined = proc.stdout + proc.stderr
            version: str | None = _extract_version(combined)
        except Exception:
            version = "unknown"
        results[name] = {"found": True, "version": version, "path": path}

    return results
