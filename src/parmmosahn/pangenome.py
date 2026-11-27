import logging
import re
import shutil
from pathlib import Path

from parmmosahn.utils import (
    ToolExecutionError,
    atomic_write_text,
    run_command,
    safe_init_output_dir,
    )

logger = logging.getLogger(__name__)


_GBK_EXTENSIONS = (".gbff", ".gbk", ".gb")


def _summarize_ppanggolin_log(log_path: Path) -> None:
    """Scan PPanGGOLiN log for notable warnings and print a brief summary."""
    if not log_path.exists():
        return
    text = log_path.read_text(encoding="utf-8", errors="replace")
    transl_except_count = text.count("transl_except")
    if transl_except_count:
        logger.info(
            f"PPanGGOLiN: {transl_except_count} transl_except warnings "
            f"(see {log_path})"
        )
    logger.info(f"Full PPanGGOLiN output saved to {log_path}")
    
    return None


def _parse_genome_count(info_output: str | None) -> int:
    """Extract the genome count from ``ppanggolin info --content`` output.

    PPanGGOLiN 2.2.x prints a ``Genomes: N`` line. We anchor on that key/value
    shape so an unrelated line that merely mentions the word (e.g. ``Genomes
    with no genes: 0``) cannot be mistaken for the count. Falls back to the
    first integer on any ``Genomes`` line for other formats. Returns 0 if no
    plausible count is found.
    """
    if not info_output:
        return 0
    # Preferred: an explicit "Genomes: N" / "Genomes = N" key/value line.
    m = re.search(r"\bGenomes\b\s*[:=]\s*([\d,]+)", info_output)
    if m:
        return int(m.group(1).replace(",", ""))
    # Fallback for older/edge formats: first integer on a "Genomes" line.
    for line in info_output.splitlines():
        if "Genomes" in line:
            digits = re.search(r"(\d[\d,]*)", line)
            if digits:
                return int(digits.group(1).replace(",", ""))
            
    return 0


def get_genbank_files(genomes_dir: Path) -> list[Path]:
    """Helper to find all genbank extensions."""
    files: list[Path] = []
    for ext in _GBK_EXTENSIONS:
        files.extend(genomes_dir.glob(f'*{ext}'))
        
    return files


def clean_genbank_files(genomes_dir: Path, work_dir: Path | None = None) -> list[Path]:
    """Return the GenBank files to feed PPanGGOLiN, cleaning curly apostrophes.

    Non-destructive: files that need cleaning are rewritten as copies in
    *work_dir* and the user's originals in *genomes_dir* are never modified.
    Clean files are returned by their original path (no copy is made). If
    *work_dir* is ``None`` no copies are written and the original paths are
    returned unchanged.
    """
    logger.info(
        f"Checking for problematic characters in GenBank files in {genomes_dir}..."
    )
    gb_files = get_genbank_files(genomes_dir)

    if not gb_files:
        logger.warning(f"No .gbff, .gbk, or .gb files found in {genomes_dir}.")
        return []

    result: list[Path] = []
    n_cleaned = 0
    for gb_path in gb_files:
        try:
            content = gb_path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            logger.error(
                f"Encoding error reading {gb_path}. Ensure it is a valid text file."
            )
            result.append(gb_path)  # still hand the original to PPanGGOLiN
            continue

        if "’" in content and work_dir is not None:
            work_dir.mkdir(parents=True, exist_ok=True)
            cleaned_path = work_dir / gb_path.name
            atomic_write_text(cleaned_path, content.replace("’", ""))
            result.append(cleaned_path)
            n_cleaned += 1
        else:
            result.append(gb_path)

    if n_cleaned:
        logger.info(
            f"Cleaned {n_cleaned} GenBank file(s) into {work_dir} "
            f"(originals left untouched)."
        )
        
    return result


def create_annotations_list(
    genomes_dir: Path, output_list: Path, gb_files: list[Path] | None = None
) -> None:
    """Generate the annotations.list file PPanGGOLiN needs, listing GenBank files.

    If *gb_files* is given those exact paths are listed (e.g. cleaned copies);
    otherwise every GenBank file under *genomes_dir* is listed.
    """
    logger.info(f"Generating annotations list at {output_list}...")
    if gb_files is None:
        gb_files = get_genbank_files(genomes_dir)

    with output_list.open("w", encoding="utf-8") as f:
        for gb in gb_files:
            f.write(f"{gb.name}\t{gb.resolve()}\n")
            
    return None


def build_pangenome(
    genomes_dir: str | Path,
    output_dir: str | Path,
    label: str,
    cluster_id: float = 0.80,
    cluster_cov: float = 0.60,
    threads: int = 1,
    tmp_dir: str | Path | None = None,
    add_small: bool = False,
    occurrence: float | None = None,
    force: bool = False,
) -> Path:
    """
    Executes the PPanGGOLiN pipeline to build a pangenome and extract the included gene families.

    Args:
        genomes_dir: Directory containing input .gbff files.
        output_dir: Directory where PPanGGOLiN will write its output.
        label: Name label for the final fasta file (e.g., 'clostridiales').
        cluster_id: Identity threshold for clustering.
        cluster_cov: Coverage threshold for clustering.
        threads: Number of CPU threads to use.
        tmp_dir: Optional temporary directory for heavy clustering operations.
        add_small: If True, export *all* gene families
        occurrence: Value in (0, 1] to override the default 2/n_genomes
             gene family occurrence threshold.
        force: If True, force overwriting any existing output

    Returns:
        Path to the exported protein-family FASTA: only included families by
        default, or the full pangenome proteome when ``add_small`` is set.
    """
    genomes_dir = Path(genomes_dir)
    output_dir = Path(output_dir)
    
    safe_init_output_dir(output_dir, force = force)

    annotations_list = output_dir / "annotations.list"
    h5_db = output_dir / "pangenome.h5"
    ppanggolin_log = output_dir / "ppanggolin.log"

    # Step 1: Pre-process inputs (non-destructive: any file needing cleaning is
    # copied into cleaned_genomes/ and the user's originals are left untouched).
    cleaned_dir = output_dir / "cleaned_genomes"
    gb_files = clean_genbank_files(genomes_dir, cleaned_dir)
    create_annotations_list(genomes_dir, annotations_list, gb_files=gb_files)

    # Step 2: Annotate (Build DB)
    logger.info("Running PPanGGOLiN annotate...")
    annotate_cmd = [
        "ppanggolin",
        "annotate",
        "-f",
        "--anno",
        str(annotations_list),
        "-o",
        str(output_dir),
        "-c",
        str(threads),
    ]
    run_command(annotate_cmd, stdout_file=ppanggolin_log)

    # Step 3: Cluster
    logger.info("Running PPanGGOLiN cluster...")
    cluster_cmd = [
        "ppanggolin",
        "cluster",
        "-f",
        "-p",
        str(h5_db),
        "--identity",
        str(cluster_id),
        "--coverage",
        str(cluster_cov),
        "-c",
        str(threads),
    ]
    if tmp_dir:
        cluster_cmd.extend(["--tmpdir", str(tmp_dir)])
    run_command(cluster_cmd, stdout_file=ppanggolin_log)

    # Step 4: Extract total genome count to calculate occurrence threshold
    logger.info("Extracting pangenome statistics...")
    info_out = run_command(["ppanggolin", "info", "-p", str(h5_db), "--content"])

    no_genomes = _parse_genome_count(info_out)
    if no_genomes < 1:
        raise ToolExecutionError(
            "Could not determine a valid genome count from `ppanggolin info "
            f"--content`. Parsed {no_genomes}. Raw output:\n{info_out}"
        )

    # Users can override this via occurrence_override if they want a different fraction.
    if occurrence is not None:
        if not 0.0 < occurrence <= 1.0:
            raise ValueError(
                f"occurrence_override must be a fraction in (0, 1], "
                f"got {occurrence}."
            )
        occurrence_threshold = round(occurrence, 5)
        logger.info(
            f"Using user-provided occurrence threshold: {occurrence_threshold} (override)."
        )
    else:
        occurrence_threshold = round(2 / no_genomes, 5)
        logger.info(
            f"Calculated occurrence threshold: {occurrence_threshold} for {no_genomes} genomes."
        )

    partitions = "all" if add_small else "softcore"

    # Step 5: Export representatives FASTA
    logger.info("Extracting representative protein families...")
    fasta_cmd = [
        "ppanggolin",
        "fasta",
        "-f",
        "-p",
        str(h5_db),
        "-o",
        str(output_dir),
        "--prot_families",
        partitions,
        "--soft_core",
        str(occurrence_threshold),
    ]
    run_command(fasta_cmd, stdout_file=ppanggolin_log)

    # Summarize notable PPanGGOLiN warnings for the user
    _summarize_ppanggolin_log(ppanggolin_log)

    # Step 6: Validate, then rename the final output file
    extracted_fasta = output_dir / f"{partitions}_protein_families.faa"
    final_fasta = output_dir / f"{label}.faa"

    if not extracted_fasta.exists():
        raise ToolExecutionError(
            f"Expected output file {extracted_fasta} was not generated."
        )

    # Guard against a degenerate representatives set (e.g. a threshold that admits no
    # family): an empty FASTA would otherwise flow silently into build_pan_model
    # and fail much later with a confusing error. Counted by streaming so a large
    # "all" proteome is not read into memory.
    n_families = 0
    with extracted_fasta.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                n_families += 1
    if n_families == 0:
        raise ToolExecutionError(
            f"PPanGGOLiN produced an empty protein-family FASTA "
            f"({extracted_fasta.name}) at occurrence threshold {occurrence_threshold}: "
            f"no families met the threshold. Try a smaller --occurrence value."
        )

    shutil.move(extracted_fasta, final_fasta)
    logger.info(
        f"Pangenome construction complete ({n_families} families)! "
        f"File saved to: {final_fasta}"
    )

    return final_fasta
