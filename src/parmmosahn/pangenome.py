import logging
import math
import shutil
import subprocess
from pathlib import Path

# Configure basic logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class ToolExecutionError(Exception):
    """Custom exception raised when an external bio-tool fails."""

    pass


def clean_genbank_files(genomes_dir: Path) -> None:
    """
    Removes problematic characters (like '’') from GenBank files.
    Replaces the bash command: find ... | xargs ... tr -d "’"
    """
    logger.info(
        f"Checking for problematic characters in GBFF files in {genomes_dir}..."
    )
    gbff_files = list(genomes_dir.glob("*.gbff"))

    if not gbff_files:
        logger.warning(f"No .gbff files found in {genomes_dir}.")
        return

    for gbff_path in gbff_files:
        try:
            content = gbff_path.read_text(encoding="utf-8")
            if "’" in content:
                logger.info(f"Cleaning problematic characters from {gbff_path.name}")
                clean_content = content.replace("’", "")
                gbff_path.write_text(clean_content, encoding="utf-8")
        except UnicodeDecodeError:
            logger.error(
                f"Encoding error reading {gbff_path}. Ensure it is a valid text file."
            )


def create_annotations_list(genomes_dir: Path, output_list: Path) -> None:
    """
    Creates a tab-separated list mapping filenames to their absolute paths.
    Replaces the bash command: paste <(dir -1 ...) <(find ...) > annotations.list
    """
    logger.info(f"Generating annotations list at {output_list}...")
    gbff_files = list(genomes_dir.glob("*.gbff"))

    with output_list.open("w", encoding="utf-8") as f:
        for gbff in gbff_files:
            # PPanGGOLiN expects: "filename \t absolute_filepath"
            f.write(f"{gbff.name}\t{gbff.resolve()}\n")


def run_command(cmd: list[str]) -> str:
    """
    Helper function to run external commands robustly.
    Captures stdout and stderr, raising an informative error on failure.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        logger.error(f"Standard Error output:\n{e.stderr}")
        raise ToolExecutionError(f"Tool failed with exit code {e.returncode}") from e


def build_pangenome(
    genomes_dir: str | Path,
    output_dir: str | Path,
    label: str,
    cluster_id: float = 0.80,
    cluster_cov: float = 0.60,
    threads: int = 1,
    tmp_dir: str | Path | None = None,
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
    for line in info_out.splitlines():
        if "Genomes" in line:
            # E.g., extracts 205 from "Genomes : 205"
            no_genomes = int("".join(filter(str.isdigit, line)))
            break

    if no_genomes == 0:
        raise ToolExecutionError(
            "Could not determine the number of genomes from PPanGGOLiN info."
        )

    # Calculate soft core threshold exactly as bash `echo "scale=5; 2/$no_genomes" | bc`
    soft_core_threshold = round(2 / no_genomes, 5)
    logger.info(
        f"Calculated soft_core threshold: {soft_core_threshold} for {no_genomes} genomes."
    )

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
        "softcore",
        "--soft_core",
        str(soft_core_threshold),
    ]
    run_command(fasta_cmd)

    # Step 6: Rename final output file
    extracted_fasta = output_dir / "softcore_protein_families.faa"
    final_fasta = output_dir / f"{label}.faa"

    if extracted_fasta.exists():
        shutil.move(extracted_fasta, final_fasta)
        logger.info(f"Pangenome construction complete! File saved to: {final_fasta}")
    else:
        raise ToolExecutionError(
            f"Expected output file {extracted_fasta} was not generated."
        )

    return final_fasta
