import contextlib
import logging
import os
import platform
import shutil
import subprocess
import tempfile
import threading
from collections import deque
from pathlib import Path
from typing import Any, Generator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class ToolExecutionError(Exception):
    pass

class FolderExistsError(FileExistsError):
    pass


def atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> Path:
    """Write *text* to *path* atomically via a temp file in the same dir.

    A crash or interrupt mid-write would otherwise leave a truncated file that
    a resumed run could mistake for complete output. Writing to a sibling temp
    file and ``os.replace``-ing it onto the target is atomic on POSIX, so the
    final path is either the old content or the new content, never a partial.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding=encoding) as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise
    return path


def safe_init_output_dir(output: Path | str, force: bool = False) -> None:
    output_path = Path(output)
    try:
        output_path.mkdir(parents=True, exist_ok=force)
    except FileExistsError:
        raise FolderExistsError(
            f'Output folder {output_path} already exists! Use --force to overwrite it.'
            )
        
        
def safe_init_output_file(output_file: Path | str, force: bool = False) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.is_file() and not force:
        raise FileExistsError(
            f'Output file {output_file} already exists! Use --force to overwrite it.'
            )

def get_system_resources() -> dict[str, Any]:
    """Best-effort inventory of host CPU, RAM, and OS for run provenance.

    Returns a dict with ``cpu_count``, ``total_ram_gb``, ``available_ram_gb``,
    ``platform``, and ``python_version``. Any field that cannot be probed is
    omitted. Memory is read from ``/proc/meminfo`` on Linux; non-Linux hosts
    will return without memory fields rather than failing.
    """
    info: dict[str, Any] = {
        "platform": platform.platform(),
        "python_version": platform.python_version(),
    }
    cpus = os.cpu_count()
    if cpus is not None:
        info["cpu_count"] = cpus

    meminfo_path = Path("/proc/meminfo")
    if meminfo_path.exists():
        try:
            kv: dict[str, int] = {}
            for line in meminfo_path.read_text().splitlines():
                k, _, rest = line.partition(":")
                tok = rest.strip().split()
                if tok and tok[0].isdigit():
                    kv[k.strip()] = int(tok[0])  # kB
            if "MemTotal" in kv:
                info["total_ram_gb"] = round(kv["MemTotal"] / 1024 / 1024, 2)
            # MemAvailable is more useful than MemFree (accounts for reclaimable cache)
            if "MemAvailable" in kv:
                info["available_ram_gb"] = round(kv["MemAvailable"] / 1024 / 1024, 2)
        except Exception:
            pass

    return info


def format_elapsed(seconds: float) -> str:
    """Format elapsed seconds into a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s}s"
    h, m = divmod(m, 60)
    return f"{h}h {m}m {s}s"


def run_command(
    cmd: list[str],
    stdout_file: Path | None = None,
    stderr_file: Path | None = None,
    timeout: float | None = None,
) -> str | None:
    """
    Runs an external command robustly with REAL-TIME logging.
    Separates stdout and stderr completely to avoid stream pollution.

    If *timeout* (seconds) is given and the process has not finished in time,
    it is killed and a ``ToolExecutionError`` is raised. ``None`` (the default)
    waits indefinitely. Pass a generous value for long external tools to guard
    against a wedged subprocess hanging the whole pipeline.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")

    captured_stdout: list[str] = []
    # Bound stderr retention: it is only used to build an error message, so a
    # tail is enough. Without this a chatty tool (× parallel strains) could
    # accumulate unbounded output in RAM. stdout is left unbounded because
    # callers consume the return value (e.g. parsing `ppanggolin info`).
    captured_stderr: deque[str] = deque(maxlen=1000)

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
                try:
                    process.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                    t_out.join()
                    t_err.join()
                    raise ToolExecutionError(
                        f"Tool timed out after {timeout}s and was killed: "
                        f"{' '.join(cmd)}"
                    )
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

        except FileNotFoundError as e:
            # Executable not on PATH — surface as the contract error every
            # caller already expects, not a raw FileNotFoundError traceback.
            raise ToolExecutionError(
                f"Tool not found: {cmd[0]!r} — is it installed and on PATH?"
            ) from e

        except OSError as e:
            raise ToolExecutionError(f"Failed to launch {cmd[0]!r}: {e}") from e


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

    atomic_write_text(
        output_tsv,
        "\n".join(f"{k}\t{v}" for k, v in mapping.items()) + "\n",
    )
    logger.debug(
        f"Gene ID mapping: {len(mapping)} RefSeq accessions from {gbff_path.name}"
    )
    return output_tsv


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


def check_bundled_data() -> dict[str, dict[str, Any]]:
    """Verify that runtime-required TSV lookup tables are reachable.

    Each entry of the form ``{found: bool, size_bytes: int | None,
    rows: int | None, path: str | None}``. Catches the failure mode where
    a wheel is installed but a ``package-data`` / ``MANIFEST.in`` entry
    is missing — the rest of the pipeline then runs until the first
    annotation step before crashing with ``FileNotFoundError``. Running
    ``parmmosahn doctor`` after install surfaces this before any
    pipeline work is done.
    """
    from importlib.resources import as_file, files

    required = (
        "mnxr_to_rhea.tsv",
        "seed_to_rhea.tsv",
        "seed_to_bigg.tsv",
        "seed_rxn_to_bigg.tsv",
    )
    results: dict[str, dict[str, Any]] = {}
    for name in required:
        info: dict[str, Any] = {
            "found": False,
            "size_bytes": None,
            "rows": None,
            "path": None,
        }
        try:
            resource = files("parmmosahn.data").joinpath(name)
            with as_file(resource) as path:
                if path.is_file():
                    info["found"] = True
                    info["size_bytes"] = path.stat().st_size
                    info["path"] = str(path)
                    # Row count is a cheap sanity check against truncation.
                    with path.open("rb") as fh:
                        info["rows"] = sum(1 for _ in fh)
        except (FileNotFoundError, ModuleNotFoundError):
            pass
        results[name] = info
    return results


def check_environment() -> dict[str, dict[str, str | bool | None]]:
    """Check for required external tools.

    Probes each tool with its ``--version`` flag (30 s timeout to accommodate
    Python-based tools like PPanGGOLiN that have slow startup import times).

    Returns:
        Nested dict of ``{tool_name: {found: bool, version: str | None, path: str | None}}``.
    """
    tools: dict[str, list[str]] = {
        "ppanggolin": ["ppanggolin", "--version"],
        "gapseq": ["gapseq", "-v"],
        "diamond": ["diamond", "--version"],
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
