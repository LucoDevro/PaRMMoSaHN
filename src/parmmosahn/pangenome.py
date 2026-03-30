import logging
import math
import re
import shutil
import subprocess
from pathlib import Path

from parmmosahn.utils import ToolExecutionError, run_command

logger = logging.getLogger(__name__)


def get_genbank_files(genomes_dir: Path) -> list[Path]:
    """Helper to find all genbank extensions."""
    files: list[Path] = []
    for ext in ["*.gbff", "*.gbk", "*.gb"]:
        files.extend(genomes_dir.glob(ext))
    return files


def clean_genbank_files(genomes_dir: Path) -> None:
    """Checks for and removes problematic characters (e.g., curly apostrophes) from GenBank files."""
    logger.info(
        f"Checking for problematic characters in GenBank files in {genomes_dir}..."
    )
    gb_files = get_genbank_files(genomes_dir)

    if not gb_files:
        logger.warning(f"No .gbff, .gbk, or .gb files found in {genomes_dir}.")
        return

    for gb_path in gb_files:
        try:
            content = gb_path.read_text(encoding="utf-8")
            if "’" in content:
                logger.info(f"Cleaning problematic characters from {gb_path.name}")
                clean_content = content.replace("’", "")
                gb_path.write_text(clean_content, encoding="utf-8")
        except UnicodeDecodeError:
            logger.error(
                f"Encoding error reading {gb_path}. Ensure it is a valid text file."
            )


def create_annotations_list(genomes_dir: Path, output_list: Path) -> None:
    """Generates the annotations.list file required by PPanGGOLiN, listing all GenBank files."""
    logger.info(f"Generating annotations list at {output_list}...")
    gb_files = get_genbank_files(genomes_dir)

    with output_list.open("w", encoding="utf-8") as f:
        for gb in gb_files:
            f.write(f"{gb.name}\t{gb.resolve()}\n")


def build_pangenome(
    genomes_dir: str | Path,
    output_dir: str | Path,
    label: str,
    cluster_id: float = 0.80,
    cluster_cov: float = 0.60,
    threads: int = 1,
    tmp_dir: str | Path | None = None,
    add_unique: bool = False,
) -> Path:
    """
    Executes the PPanGGOLiN pipeline to build a pangenome and extract the soft-core.

    Args:
        genomes_dir: Directory containing input .gbff files.
        output_dir: Directory where PPanGGOLiN will write its output.
        label: Name label for the final fasta file (e.g., 'clostridiales').
        cluster_id: Identity threshold for clustering.
        cluster_cov: Coverage threshold for clustering.
        threads: Number of CPU threads to use.
        tmp_dir: Optional temporary directory for heavy clustering operations.

    Returns:
        Path to the final soft-core protein families FASTA file.
    """
    genomes_dir = Path(genomes_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    annotations_list = output_dir / "annotations.list"
    h5_db = output_dir / "pangenome.h5"

    # Step 1: Pre-process inputs
    clean_genbank_files(genomes_dir)
    create_annotations_list(genomes_dir, annotations_list)

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
    run_command(annotate_cmd)

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
    run_command(cluster_cmd)

    # Step 4: Extract total genome count to calculate soft-core threshold
    logger.info("Extracting pangenome statistics...")
    info_out = run_command(["ppanggolin", "info", "-p", str(h5_db), "--content"])

    no_genomes = 0
    if info_out:
        for line in info_out.splitlines():
            if "Genomes" in line:
                match = re.search(r"\d+", line)
                if match:
                    no_genomes = int(match.group())
                break

    if no_genomes == 0:
        raise ToolExecutionError(
            "Could not determine the number of genomes from PPanGGOLiN info."
        )

    # Calculate soft core threshold exactly as bash `echo "scale=5; 2/$no_genomes" | bc`
    # TODO check if this is biologically relevant
    soft_core_threshold = round(2 / no_genomes, 5)
    logger.info(
        f"Calculated soft_core threshold: {soft_core_threshold} for {no_genomes} genomes."
    )

    # Option to include unique/cloud genes
    partitions = "all" if add_unique else "softcore"

    # Step 5: Export Soft-core FASTA
    logger.info("Extracting soft-core protein families...")
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
        str(soft_core_threshold),
    ]
    run_command(fasta_cmd)

    # Step 6: Rename final output file
    extracted_fasta = output_dir / f"{partitions}_protein_families.faa"
    final_fasta = output_dir / f"{label}.faa"

    if extracted_fasta.exists():
        shutil.move(extracted_fasta, final_fasta)
        logger.info(f"Pangenome construction complete! File saved to: {final_fasta}")
    else:
        raise ToolExecutionError(
            f"Expected output file {extracted_fasta} was not generated."
        )

    return final_fasta
