import logging
import re
import shutil
from pathlib import Path

import pandas as pd

from parmmosahn.utils import ToolExecutionError, extract_proteome_from_gbff, run_command

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
        "-M",
        "prot",  # Explicitly tell gapseq these are proteins
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
    Gracefully accepts both .gbff (extracts CDS) and .faa (direct mapping).
    """
    genomes_dir = Path(genomes_dir).resolve()
    pangenome_fasta = Path(pangenome_fasta).resolve()
    panmodel_rxn_table = Path(panmodel_rxn_table).resolve()
    medium_file = Path(medium_file).resolve()
    output_dir = Path(output_dir).resolve()

    # HPC Readiness: Create an isolated directory for processed proteomes
    matches_dir = output_dir / "matches"
    models_dir = output_dir / "models"
    proteomes_dir = output_dir / "proteomes"

    for d in [matches_dir, models_dir, proteomes_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Collect both file types
    faa_files = list(genomes_dir.glob("*.faa"))
    gbff_files = list(genomes_dir.glob("*.gbff"))

    if not faa_files and not gbff_files:
        raise FileNotFoundError(f"No .faa or .gbff files found in {genomes_dir}")

    processed_genomes = []

    # 1. Prepare Genomes (Extract GBFFs and isolate FAAs)
    logger.info("Preparing input genomes for DIAMOND mapping...")

    for gbff in gbff_files:
        out_faa = proteomes_dir / f"{gbff.stem}.faa"
        extract_proteome_from_gbff(gbff, out_faa)
        clean_fasta_file(out_faa)  # Clean the safe, isolated file
        processed_genomes.append(out_faa)

    for faa in faa_files:
        out_faa = proteomes_dir / faa.name
        if not out_faa.exists():
            shutil.copy(faa, out_faa)  # Protect read-only original inputs
        clean_fasta_file(out_faa)
        processed_genomes.append(out_faa)

    # Step 2: Build Diamond DB
    logger.info("Building Diamond database from pangenome...")
    db_path = output_dir / pangenome_fasta.stem
    run_command(["diamond", "makedb", "--in", str(pangenome_fasta), "-d", str(db_path)])
    db_file = str(db_path) + ".dmnd"

    # Step 3: Iterate through strains
    for genome in processed_genomes:
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

        # 3c. Find Transporters (Using the protein FASTA)
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
        final_xml = strain_dir / f"{strain_name}.xml"

        if xml_file.exists():
            xml_content = xml_file.read_text(encoding="utf-8")

            # 1. Fix the exchange reaction tags
            if "R_EX_" in xml_content:
                xml_content = xml_content.replace("R_EX_", "EX_")

            # 2. Fix the internal SBML Model ID (Sanitize for SBML SId rules)
            # Replace any character that is NOT a letter, number, or underscore with an underscore
            safe_model_id = re.sub(r"[^a-zA-Z0-9_]", "_", strain_name)

            # SBML IDs also cannot start with a number. Catch that edge case just in case.
            if safe_model_id and safe_model_id[0].isdigit():
                safe_model_id = f"M_{safe_model_id}"

            xml_content = xml_content.replace(
                'id="filtered_Reactions_tbl"', f'id="{safe_model_id}"'
            )
            xml_content = xml_content.replace(
                'id="filtered-Reactions.tbl"', f'id="{safe_model_id}"'
            )

            # Write to the new properly named file and delete the old gapseq one
            final_xml.write_text(xml_content, encoding="utf-8")
            xml_file.unlink()
            logger.debug(
                f"Saved cleaned XML as {final_xml.name} with internal ID {safe_model_id}"
            )

        # This prevents memote and the pipeline from picking up duplicate/unfinished models
        draft_xml = strain_dir / "filtered-Reactions.tbl-draft.xml"
        if draft_xml.exists():
            draft_xml.unlink()
            logger.debug(f"Removed intermediate draft model {draft_xml.name}")

    logger.info(f"All strain models successfully generated in {models_dir}")
    return models_dir
