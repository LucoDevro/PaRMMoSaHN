import logging
from pathlib import Path

import pandas as pd

from parmmosahn.utils import ToolExecutionError, run_command

logger = logging.getLogger(__name__)


def clean_fasta_file(fasta_path: Path) -> None:
    """
    Removes problematic characters (like '’') from FASTA files to prevent gapseq errors.
    Replaces: cat $pangenome | tr -d "’" > tmp && mv tmp $pangenome
    """
    logger.info(f"Cleaning problematic characters from {fasta_path.name}...")
    try:
        content = fasta_path.read_text(encoding="utf-8")
        if "’" in content:
            clean_content = content.replace("’", "")
            fasta_path.write_text(clean_content, encoding="utf-8")
            logger.debug(f"Cleaned {fasta_path.name}.")
    except UnicodeDecodeError:
        logger.error(
            f"Encoding error reading {fasta_path}. Ensure it is a valid text file."
        )
        raise


def build_pan_model(
    pangenome_fasta: str | Path,
    output_dir: str | Path,
    blast_sens: int = 150,
    threads: int = 14,
) -> Path:
    """
    Executes 'gapseq find' on the pangenome FASTA to construct a pan-metabolic model.

    Args:
        pangenome_fasta: Path to the soft-core FASTA from PPanGGOLiN.
        output_dir: Directory to save the gapseq output and logs.
        blast_sens: Blast sensitivity parameter for gapseq (default: 150).
        threads: Number of CPU cores to use.

    Returns:
        Path to the output directory containing the pan-model results.
    """
    fasta_path = Path(pangenome_fasta).resolve()
    out_dir = Path(output_dir).resolve()

    if not fasta_path.exists():
        raise FileNotFoundError(f"Pangenome FASTA not found at: {fasta_path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Prepare input
    clean_fasta_file(fasta_path)

    # 2. Setup log files
    log_out = out_dir / "log.log"
    log_err = out_dir / "log.err"

    # 3. Construct gapseq command
    cmd = [
        "gapseq",
        "find",
        "-p",
        "all",
        "-v",
        "1",
        "-y",
        "-j",
        "-O",
        "-b",
        str(blast_sens),
        "-f",
        str(out_dir),  # Gapseq output directory
        "-K",
        str(threads),
        str(fasta_path),
    ]

    logger.info(
        "Running gapseq find. This may take a while depending on the pangenome size..."
    )

    # 4. Execute (outputs are piped directly to the log files, keeping RAM usage low)
    run_command(cmd, stdout_file=log_out, stderr_file=log_err)

    logger.info(f"Pan-model generation complete! Results saved in {out_dir}")

    return out_dir


def filter_blast_table(
    pan_reaction_table_file: Path, strain_blast_hits_file: Path, output_file: Path
) -> None:
    """
    Filters the panmodel reaction table based on DIAMOND blast hits.
    Native implementation of the provided 'filter_blast_table.py'.
    """

    def format_number(num):
        return "{:.0f}".format(num) if num.is_integer() else str(num)

    logger.debug(f"Filtering reactions based on {strain_blast_hits_file.name}")

    # Read diamond output
    strain_blast_hits = pd.read_table(
        strain_blast_hits_file,
        header=None,
        usecols=[0, 1, 11],
        names=["Query", "Target", "Repr. bitscore"],
    )

    # Read pan-model reactions
    pan_reaction_table = pd.read_table(pan_reaction_table_file, comment="#")

    # Merge and filter
    strain_reaction_table = pd.merge(
        strain_blast_hits,
        pan_reaction_table,
        how="inner",
        left_on="Target",
        right_on="stitle",
    )

    strain_reaction_table.drop_duplicates(
        inplace=True, ignore_index=True, subset=["Query", "Target", "rxn", "pathway"]
    )
    strain_reaction_table.sort_values(
        by="Repr. bitscore", ascending=False, inplace=True
    )
    strain_reaction_table.drop_duplicates(
        subset=["Query", "pathway"], keep="first", inplace=True
    )

    strain_reaction_table.sort_index(inplace=True)
    strain_reaction_table["stitle"] = strain_reaction_table["Query"]
    strain_reaction_table.drop(
        columns=["Query", "Target", "Repr. bitscore"], inplace=True
    )

    # Save output
    strain_reaction_table.to_csv(
        output_file,
        sep="\t",
        index=False,
        na_rep="NA",
        float_format=format_number,
    )


def derive_strain_models(
    genomes_dir: str | Path,
    pangenome_fasta: str | Path,
    panmodel_rxn_table: str | Path,
    output_dir: str | Path,
    medium_file: str | Path,
    diamond_bits: int = 150,
    gapseq_bits: int = 150,
    biomass: str = "pos",
    threads: int = 14,
) -> Path:
    """
    Uses diamond and gapseq to project the pan-model onto individual strains.
    Replaces 03-derive_strain_models.ipynb logic.
    """
    genomes_dir = Path(genomes_dir).resolve()
    pangenome_fasta = Path(pangenome_fasta).resolve()
    panmodel_rxn_table = Path(panmodel_rxn_table).resolve()
    medium_file = Path(medium_file).resolve()
    output_dir = Path(output_dir).resolve()

    matches_dir = output_dir / "matches"
    models_dir = output_dir / "models"
    matches_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)

    genomes = list(genomes_dir.glob("*.faa"))
    if not genomes:
        raise FileNotFoundError(f"No .faa genome files found in {genomes_dir}")

    # Step 1: Clean FASTA files
    for genome in genomes:
        clean_fasta_file(genome)  # Using the function we built in Step 2

    # Step 2: Build Diamond DB
    logger.info("Building Diamond database from pangenome...")
    db_path = output_dir / pangenome_fasta.stem
    run_command(["diamond", "makedb", "--in", str(pangenome_fasta), "-d", str(db_path)])

    # We must append .dmnd for the blastp command
    db_file = str(db_path) + ".dmnd"

    # Step 3: Iterate through strains to map and generate models
    for genome in genomes:
        strain_name = genome.stem
        logger.info(f"--- Processing strain: {strain_name} ---")

        strain_dir = models_dir / strain_name
        strain_dir.mkdir(exist_ok=True)

        log_out = strain_dir / "log.log"
        log_err = strain_dir / "log.err"

        # 3a. Diamond Blastp
        match_file = matches_dir / f"{strain_name}.tsv"
        run_command(
            [
                "diamond",
                "blastp",
                "-d",
                db_file,
                "-q",
                str(genome),
                "-o",
                str(match_file),
                "--min-score",
                str(diamond_bits),
                "-p",
                str(threads),
            ]
        )

        # 3b. Filter Panmodel table
        filtered_rxn_tbl = strain_dir / "filtered-Reactions.tbl"
        filter_blast_table(panmodel_rxn_table, match_file, filtered_rxn_tbl)

        # 3c. Find Transporters
        run_command(
            [
                "gapseq",
                "find-transport",
                "-K",
                str(threads),
                "-b",
                str(gapseq_bits),
                "-f",
                str(strain_dir),
                str(genome),
            ],
            stdout_file=log_out,
            stderr_file=log_err,
        )

        # 3d. Draft Model
        transporter_tbl = strain_dir / f"{strain_name}-Transporter.tbl"
        run_command(
            [
                "gapseq",
                "draft",
                "-r",
                str(filtered_rxn_tbl),
                "-t",
                str(transporter_tbl),
                "-b",
                biomass,
                "-u",
                str(gapseq_bits),
                "-f",
                str(strain_dir),
            ],
            stdout_file=log_out,
            stderr_file=log_err,
        )

        # 3e. Gapfill Model
        draft_rds = strain_dir / "filtered-Reactions.tbl-draft.RDS"
        weights_rds = strain_dir / "filtered-Reactions.tbl-rxnWeights.RDS"
        genes_rds = strain_dir / "filtered-Reactions.tbl-rxnXgenes.RDS"

        run_command(
            [
                "gapseq",
                "fill",
                "-q",
                "-m",
                str(draft_rds),
                "-n",
                str(medium_file),
                "-c",
                str(weights_rds),
                "-g",
                str(genes_rds),
                "-f",
                str(strain_dir),
            ],
            stdout_file=log_out,
            stderr_file=log_err,
        )

        # 3f. Clean XML tags (Replaces the bash `sed` command)
        xml_file = strain_dir / "filtered-Reactions.tbl.xml"
        if xml_file.exists():
            xml_content = xml_file.read_text(encoding="utf-8")
            if "R_EX_" in xml_content:
                xml_file.write_text(
                    xml_content.replace("R_EX_", "EX_"), encoding="utf-8"
                )
                logger.debug(f"Fixed 'R_EX_' tags in {xml_file.name}")

    logger.info(f"All strain models successfully generated in {models_dir}")
    return models_dir
