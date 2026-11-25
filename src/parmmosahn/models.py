import hashlib
import logging
import re
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

from parmmosahn.utils import (
    ToolExecutionError,
    atomic_write_text,
    extract_gene_id_mapping,
    extract_proteome_from_gbff,
    run_command,
)

logger = logging.getLogger(__name__)


def clean_fasta_file(fasta_path: Path, output_path: Path | None = None) -> Path:
    """
    Clean a FASTA file by:
    1. Truncating headers at the first space (descriptions break gapseq's regex).
    2. Stripping trailing ``*`` stop codons that PPanGGOLiN adds to representative
       sequences (these make R/Biostrings report thousands of "invalid one-letter
       sequence codes" during gapseq BLAST).
    3. Removing the U+2019 right single quote left by some upstream tools.

    If *output_path* is given the cleaned FASTA is written there and *fasta_path*
    is left untouched; otherwise the file is cleaned in place. The write is atomic,
    so an interrupted run never leaves a half-cleaned file. Returns the path that
    now holds the cleaned FASTA.
    """
    dest = Path(output_path) if output_path is not None else fasta_path
    logger.debug(f"Cleaning formatting for {fasta_path.name}...")
    try:
        lines = fasta_path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        logger.error(
            f"Encoding error reading {fasta_path}. Ensure it is a valid text file."
        )
        raise

    cleaned_lines = []
    for line in lines:
        if line.startswith(">"):
            # Keep only the pure ID (e.g. >AEP89073.1)
            line = line.split()[0]
        else:
            line = line.rstrip("*")

        if "’" in line:
            line = line.replace("’", "")

        cleaned_lines.append(line)

    atomic_write_text(dest, "\n".join(cleaned_lines) + "\n")
    logger.debug(f"Cleaned {fasta_path.name}.")
    return dest


def build_pan_model(
    pangenome_fasta: str | Path,
    output_dir: str | Path,
    blast_sens: int = 150,
    threads: int = 14,
    force: bool = False,
) -> Path:
    """
    Executes 'gapseq find' on the pangenome FASTA to construct a pan-metabolic model.

    Args:
        pangenome_fasta: Path to the representatives FASTA from PPanGGOLiN.
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
        
    try:
        out_dir.mkdir(parents=True, exist_ok=force)
    except FileExistsError:
        raise ToolExecutionError(
            "Panmodel results folder already exists. Use --force to overwrite it."
            )

    # 1. Prepare input: clean into a derived copy so the source representatives FASTA
    #    is never mutated in place — a re-run or any external consumer must still
    #    see the original headers.
    cleaned_fasta = clean_fasta_file(
        fasta_path, output_path=out_dir / f"{fasta_path.stem}.faa"
    )

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
        str(cleaned_fasta),
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
    Filters the panmodel reaction table based on DIAMOND blast hits using Polars LazyFrames
    """
    logger.debug(f"Filtering reactions based on {strain_blast_hits_file.name}")

    # 1. Setup LazyFrame for Blast Hits
    blast_cols = [
        "Query",
        "Target",
        "drop1",
        "drop2",
        "drop3",
        "drop4",
        "drop5",
        "drop6",
        "drop7",
        "drop8",
        "drop9",
        "Repr. bitscore",
    ]

    # FIX: Added quote_char=None to prevent stray quotes from breaking chunk parsing
    lf_blast = (
        pl.scan_csv(
            strain_blast_hits_file,
            separator="\t",
            has_header=False,
            new_columns=blast_cols,
            quote_char=None,
        )
        .select(["Query", "Target", "Repr. bitscore"])
        .unique(subset=["Query"], keep="first")
    )

    # 2. Setup LazyFrame for Pan-model Reactions
    # FIX: Added quote_char=None, truncate_ragged_lines, and increased infer_schema_length
    df_pan = pl.read_csv(
        pan_reaction_table_file,
        separator="\t",
        comment_prefix="#",
        quote_char=None,
        truncate_ragged_lines=True,
        infer_schema_length=10000,
    )
    lf_pan = df_pan.lazy()

    # 3. Perform the joins and complex filtering lazily
    lf_filtered = (
        lf_blast.join(lf_pan, left_on="Target", right_on="stitle", how="inner")
        .unique(subset=["Query", "Target", "rxn", "pathway"], keep="first")
        .sort("Repr. bitscore", descending=True)
        .unique(subset=["Query", "pathway", "dbhit"], keep="first")
        .sort("Query")
        .with_columns(pl.col("Query").alias("stitle"))
        .drop(["Query", "Target", "Repr. bitscore"])
    )

    # 4. Trigger execution graph and stream to disk
    # FIX: Added quote_style="never" to ensure safe output for gapseq
    lf_filtered.collect().write_csv(
        output_file, separator="\t", null_value="NA", quote_style="never"
    )


def _process_single_strain(
    genome: Path,
    total_strains: int,
    models_workdir: Path,
    matches_dir: Path,
    db_file: Path,
    panmodel_rxn_table: Path,
    medium_file: Path,
    strain_bits: int,
    pan_bits: int,
    biomass: str,
    threads_per_strain: int,
    gapfill: bool,
    counter: "list[int]",
    counter_lock: threading.Lock,
    model_id_map: dict[str, str],
) -> tuple[str, bool]:
    """Process a single strain end-to-end. Returns (strain_name, skipped)."""
    strain_name = genome.stem
    strain_dir = models_workdir / strain_name
    strain_dir.mkdir(exist_ok=True)
    done_marker = strain_dir / ".done"

    with counter_lock:
        counter[0] += 1
        display_idx = counter[0]

    if done_marker.exists():
        logger.info(
            f"[{display_idx}/{total_strains}] Skipping {strain_name} (already completed)"
        )
        return strain_name, True

    logger.info(
        f"[{display_idx}/{total_strains}] --- Processing strain: {strain_name} ---"
    )

    log_out = strain_dir / "log.log"
    log_err = strain_dir / "log.err"

    # 3a. Diamond Blastp
    match_file = matches_dir / f"{strain_name}.tsv"
    run_command(
        [
            "diamond",
            "blastp",
            "-d",
            str(db_file),
            "-q",
            str(genome),
            "-o",
            str(match_file),
            "--min-score",
            str(strain_bits),
            "-p",
            str(threads_per_strain),
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
            str(threads_per_strain),
            "-b",
            str(pan_bits),
            "-M",
            "prot",
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
            str(pan_bits),
            "-f",
            str(strain_dir),
        ],
        stdout_file=log_out,
        stderr_file=log_err,
    )

    # 3e. Gapfill Model (or pass through the draft model)
    draft_xml = strain_dir / "filtered-Reactions.tbl-draft.xml"
    draft_rds = strain_dir / "filtered-Reactions.tbl-draft.RDS"
    weights_rds = strain_dir / "filtered-Reactions.tbl-rxnWeights.RDS"
    genes_rds = strain_dir / "filtered-Reactions.tbl-rxnXgenes.RDS"
    
    xml_file = strain_dir / "filtered-Reactions.tbl.xml"
    if gapfill:
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
    
    else:
        shutil.copy(draft_xml, xml_file)
        
    # Remove intermediate model
    if draft_xml.exists():
        draft_xml.unlink()
        logger.debug(f"Removed intermediate draft model {draft_xml.name}")

    # 3f. Clean XML
    final_xml = strain_dir / f"{strain_name}.xml"

    if xml_file.exists():
        xml_content = xml_file.read_text(encoding="utf-8")
        xml_content = xml_content.replace("R_EX_", "EX_")
        # Collision-free id assigned up front in derive_strain_models.
        safe_model_id = model_id_map[strain_name]
        xml_content = re.sub(
            r'((?:id|metaid)=")filtered[-_]Reactions[._]?tbl"',
            r"\g<1>" + safe_model_id + '"',
            xml_content,
        )
        xml_content = re.sub(
            r'(name=")filtered[-.]Reactions[._]?tbl"',
            r"\g<1>" + strain_name + '"',
            xml_content,
        )
        final_xml.write_text(xml_content, encoding="utf-8")
        xml_file.unlink()
        logger.debug(
            f"Saved cleaned XML as {final_xml.name} with internal ID {safe_model_id}"
        )

    # Only checkpoint once a real model is on disk. gapseq `fill` can exit 0
    # without producing an XML on degenerate input; a premature .done would make
    # a resumed run skip the strain forever with no model and silently shrink
    # the cohort. Raising here lets the strain be reported as failed and retried.
    if not final_xml.exists():
        raise ToolExecutionError(
            f"gapseq produced no model for {strain_name} "
            f"(expected {final_xml.name}); strain not checkpointed."
        )

    done_marker.write_text(
        datetime.now(timezone.utc).isoformat(timespec="seconds"),
        encoding="utf-8",
    )
    return strain_name, False


def derive_strain_models(
    genomes_dir: str | Path,
    pangenome_fasta: str | Path,
    panmodel_rxn_table: str | Path,
    output_dir: str | Path,
    medium_file: str | Path,
    gapfill: bool,
    strain_bits: int = 150,
    pan_bits: int = 150,
    biomass: str = "pos",
    threads: int = 14,
    parallel_strains: int | None = None,
    force: bool = False,
) -> Path:
    """
    Uses diamond and gapseq to project the pan-model onto individual strains.
    Gracefully accepts both .gbff (extracts CDS) and .faa (direct mapping).

    Args:
        parallel_strains: Number of strains to process concurrently. Each strain
            receives threads // parallel_strains CPU threads. Default 1 (sequential).
            Recommended: 2–4 on 8–16 GB machines (each strain peaks at ~1.5 GB RAM).
    """
    genomes_dir = Path(genomes_dir).resolve()
    pangenome_fasta = Path(pangenome_fasta).resolve()
    panmodel_rxn_table = Path(panmodel_rxn_table).resolve()
    medium_file = Path(medium_file).resolve()
    output_dir = Path(output_dir).resolve()

    # Create an isolated directory for processed proteomes
    matches_dir = output_dir / "matches"
    models_dir = output_dir / "models"
    models_workdir = output_dir / "models_workdir"
    proteomes_dir = output_dir / "proteomes"
    
    try:
        for d in [matches_dir, models_workdir, proteomes_dir, models_dir]:
            d.mkdir(parents=True, exist_ok=force)
    except FileExistsError:
        raise ToolExecutionError(
            "Projection results folder already exists. Use --force to overwrite it."
            )

    # Collect both file types
    faa_files = list(genomes_dir.glob("*.faa"))
    gbff_files = list(genomes_dir.glob("*.gbff"))

    if not faa_files and not gbff_files:
        raise FileNotFoundError(f"No .faa or .gbff files found in {genomes_dir}")

    # 1. Prepare Genomes (Extract GBFFs and isolate FAAs)
    logger.info("Preparing input genomes for DIAMOND mapping...")

    processed_genomes_map: dict[str, Path] = {}
    prep_failures: list[str] = []
    for gbff in gbff_files:
        try:
            out_faa = proteomes_dir / f"{gbff.stem}.faa"
            extract_proteome_from_gbff(gbff, out_faa)
            clean_fasta_file(out_faa)
            # Build Bakta protein_id → NCBI RefSeq mapping for annotation enrichment
            mapping_tsv = proteomes_dir / f"{gbff.stem}.gene_mapping.tsv"
            if not mapping_tsv.exists():
                extract_gene_id_mapping(gbff, mapping_tsv)
            processed_genomes_map[gbff.stem] = out_faa
        except Exception as exc:  # a single bad genome must not abort prep
            prep_failures.append(gbff.name)
            logger.error(f"Skipping {gbff.name}: genome preparation failed: {exc}")

    for faa in faa_files:
        if faa.stem in processed_genomes_map:
            logger.warning(
                f"Ignoring {faa.name} as {faa.stem}.gbff was already processed."
            )
            continue
        try:
            out_faa = proteomes_dir / faa.name
            if not out_faa.exists():
                shutil.copy(faa, out_faa)
            clean_fasta_file(out_faa)
            processed_genomes_map[faa.stem] = out_faa
        except Exception as exc:  # a single bad genome must not abort prep
            prep_failures.append(faa.name)
            logger.error(f"Skipping {faa.name}: genome preparation failed: {exc}")

    if prep_failures:
        logger.warning(
            f"Genome preparation failed for {len(prep_failures)} input(s): "
            f"{', '.join(prep_failures)}"
        )

    processed_genomes = list(processed_genomes_map.values())
    if not processed_genomes:
        raise RuntimeError(
            "No genomes could be prepared for strain projection "
            f"({len(gbff_files) + len(faa_files)} input(s), all failed)."
        )

    # Pre-assign collision-free SBML model ids. Sanitising non-alphanumerics to
    # "_" could map two distinct strain names to the same id; resolve any clash
    # deterministically with a short hash suffix (no change for the common case).
    model_id_map: dict[str, str] = {}
    used_ids: set[str] = set()
    for genome in processed_genomes:
        strain_name = genome.stem
        safe = re.sub(r"[^a-zA-Z0-9_]", "_", strain_name)
        if safe and safe[0].isdigit():
            safe = f"M_{safe}"
        if safe in used_ids:
            suffix = hashlib.sha1(strain_name.encode()).hexdigest()[:6]
            logger.warning(
                f"Model-id collision for {strain_name!r} (sanitised to {safe!r}); "
                f"disambiguating with suffix _{suffix}."
            )
            safe = f"{safe}_{suffix}"
        used_ids.add(safe)
        model_id_map[strain_name] = safe

    # Step 2: Build Diamond DB
    logger.info("Building Diamond database from pangenome...")
    db_path = output_dir / pangenome_fasta.stem
    run_command(["diamond", "makedb", "--in", str(pangenome_fasta), "-d", str(db_path)])
    db_file = db_path.with_suffix(".dmnd")

    # Step 3: Process strains (parallel or sequential)
    total_strains = len(processed_genomes)
    # Thread budget: `parallel_strains` workers each run gapseq/DIAMOND with
    # `threads // parallel_strains` threads, so total CPU use stays ~`threads`.
    # When unspecified, default to threads//4 workers (~4 threads each) — a
    # balance that keeps each worker usefully parallel while bounding peak RAM
    # (each strain peaks at ~1.5 GB).
    if parallel_strains is None:
        parallel_strains = max(1, threads // 4)
    threads_per_strain = max(1, threads // parallel_strains)

    if parallel_strains > 1:
        logger.info(
            f"Processing {total_strains} strains with {parallel_strains} parallel workers "
            f"({threads_per_strain} threads each)"
        )
    else:
        logger.info(
            f"Processing {total_strains} strains sequentially ({threads} threads each)"
        )

    counter: list[int] = [0]
    counter_lock = threading.Lock()
    skipped = 0
    failed = 0

    with ThreadPoolExecutor(max_workers=parallel_strains) as pool:
        futures = {
            pool.submit(
                _process_single_strain,
                genome,
                total_strains,
                models_workdir,
                matches_dir,
                db_file,
                panmodel_rxn_table,
                medium_file,
                strain_bits,
                pan_bits,
                biomass,
                threads_per_strain,
                gapfill,
                counter,
                counter_lock,
                model_id_map,
            ): genome
            for genome in processed_genomes
        }
        try:
            for future in as_completed(futures):
                try:
                    _, was_skipped = future.result()
                    if was_skipped:
                        skipped += 1
                except Exception as exc:
                    genome = futures[future]
                    failed += 1
                    logger.error(f"Strain {genome.stem} failed: {exc}")
        except KeyboardInterrupt:
            # On Ctrl-C, cancel queued (not-yet-started) strain jobs. Strains
            # already running cannot be interrupted here — their gapseq/DIAMOND
            # subprocesses must finish before the ThreadPoolExecutor drains — so
            # the prompt returns once the in-flight workers complete.
            logger.warning("Interrupted — cancelling queued strain jobs.")
            for f in futures:
                f.cancel()
            raise
            
    # Step 4: Gather all generated models into one output folder
    for strain_dir in models_workdir.iterdir():
        strain = strain_dir.name
        shutil.copy(
            str(models_workdir / strain / strain) + '.xml',
            str(models_dir / strain) + '.xml'
            )

    if skipped:
        logger.info(
            f"Resumed: skipped {skipped}/{total_strains} previously completed strains."
        )
    generated = total_strains - failed
    if failed:
        logger.warning(
            f"Strain projection finished with {failed}/{total_strains} failure(s); "
            f"{generated} model(s) available in {models_dir}"
        )
    else:
        logger.info(
            f"All {generated} strain models successfully generated in {models_dir}"
        )
    return models_dir
