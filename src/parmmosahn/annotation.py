"""Offline annotation enrichment for FAIR-compliant metabolic models.

Recovers database cross-references that gapseq and ModelPolisher don't add:
- Genes: NCBI protein accessions from mangled SId-safe gene IDs
- Reactions: BRENDA from EC codes (1:1 mapping), RHEA from SEED/MetaNetX via bundled TSV
"""

import csv
import logging
import re
from importlib.resources import files
from pathlib import Path

import cobra

logger = logging.getLogger(__name__)

# identifiers.org URI prefixes used in SBML annotation RDF
_URI_PREFIX = {
    "ncbiprotein": "https://identifiers.org/ncbiprotein:",
    "brenda": "https://identifiers.org/brenda:",
    "rhea": "https://identifiers.org/rhea:",
}

# Pattern for gapseq-mangled gene IDs:
#   {accession}_{version}_{start}_{end}  e.g. AEP93129_1_1_308
# We recover {accession}.{version} as the NCBI protein accession.
# Handles underscore-containing prefixes like NP_391935_1_1_308 → NP_391935.1
_GENE_ID_RE = re.compile(r"^(?P<acc>[A-Za-z]{2,3}_?\d+)_(?P<ver>\d+)_\d+_\d+$")


def _load_tsv_mapping(tsv_name: str) -> dict[str, str]:
    """Load a two-column TSV from the bundled data directory into a dict."""
    resource = files("parmmosahn.data").joinpath(tsv_name)
    mapping: dict[str, str] = {}
    with resource.open("r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 2:
                mapping[row[0]] = row[1]
    return mapping


def _get_annotation(entity, key: str) -> list[str]:
    """Get annotation values for a given key from a cobra entity."""
    annot = getattr(entity, "annotation", {})
    if not isinstance(annot, dict):
        return []
    val = annot.get(key)
    if val is None:
        return []
    if isinstance(val, str):
        return [val]
    if isinstance(val, list):
        return val
    return []


def _add_annotation(entity, key: str, value: str) -> bool:
    """Add a single annotation value if not already present. Returns True if added."""
    annot = getattr(entity, "annotation", None)
    if annot is None:
        entity.annotation = {}
        annot = entity.annotation

    existing = annot.get(key)
    if existing is None:
        annot[key] = value
        return True
    if isinstance(existing, str):
        if existing == value:
            return False
        annot[key] = [existing, value]
        return True
    if isinstance(existing, list):
        if value in existing:
            return False
        existing.append(value)
        return True
    return False


def annotate_genes(model: cobra.Model) -> int:
    """Recover NCBI protein accessions from gapseq-mangled gene IDs.

    Gene IDs follow the pattern {accession}_{version}_{start}_{end},
    e.g. AEP93129_1_1_308 -> AEP93129.1 (NCBI protein accession).

    Returns the number of genes annotated.
    """
    count = 0
    for gene in model.genes:
        if _get_annotation(gene, "ncbiprotein"):
            continue
        m = _GENE_ID_RE.match(gene.id)
        if m:
            accession = f"{m.group('acc')}.{m.group('ver')}"
            if _add_annotation(gene, "ncbiprotein", accession):
                count += 1
    return count


def annotate_reactions(
    model: cobra.Model,
    seed_rhea: dict[str, str] | None = None,
    mnxr_rhea: dict[str, str] | None = None,
) -> dict[str, int]:
    """Enrich reaction annotations with BRENDA (from EC) and RHEA (from SEED/MetaNetX).

    Parameters
    ----------
    model : cobra.Model
        The model to annotate in-place.
    seed_rhea : dict, optional
        Mapping of SEED reaction IDs to RHEA IDs. Loaded from bundled data if None.
    mnxr_rhea : dict, optional
        Mapping of MetaNetX reaction IDs to RHEA IDs. Loaded from bundled data if None.

    Returns
    -------
    dict with keys 'brenda' and 'rhea', each mapping to the count of annotations added.
    """
    if seed_rhea is None:
        seed_rhea = _load_tsv_mapping("seed_to_rhea.tsv")
    if mnxr_rhea is None:
        mnxr_rhea = _load_tsv_mapping("mnxr_to_rhea.tsv")

    counts = {"brenda": 0, "rhea": 0}

    for rxn in model.reactions:
        # BRENDA: direct 1:1 from EC code
        if not _get_annotation(rxn, "brenda"):
            ec_codes = _get_annotation(rxn, "ec-code")
            for ec in ec_codes:
                if _add_annotation(rxn, "brenda", ec):
                    counts["brenda"] += 1
                    break  # one BRENDA annotation per reaction

        # RHEA: via SEED reaction ID or MetaNetX reaction ID
        if not _get_annotation(rxn, "rhea"):
            # Try seed.reaction first
            seed_ids = _get_annotation(rxn, "seed.reaction")
            for sid in seed_ids:
                if sid in seed_rhea:
                    if _add_annotation(rxn, "rhea", seed_rhea[sid]):
                        counts["rhea"] += 1
                        break

            # If still no RHEA, try metanetx.reaction
            if not _get_annotation(rxn, "rhea"):
                mnx_ids = _get_annotation(rxn, "metanetx.reaction")
                for mid in mnx_ids:
                    if mid in mnxr_rhea:
                        if _add_annotation(rxn, "rhea", mnxr_rhea[mid]):
                            counts["rhea"] += 1
                            break

    return counts


def enrich_model_annotations(
    model_path: Path,
    output_path: Path,
) -> dict[str, int]:
    """Load an SBML model, enrich annotations, and write back.

    Parameters
    ----------
    model_path : Path
        Input SBML model file.
    output_path : Path
        Output SBML model file (can be same as input for in-place).

    Returns
    -------
    dict summarizing counts: 'genes_ncbiprotein', 'reactions_brenda', 'reactions_rhea'.
    """
    model = cobra.io.read_sbml_model(str(model_path))

    # Load bundled mappings once
    seed_rhea = _load_tsv_mapping("seed_to_rhea.tsv")
    mnxr_rhea = _load_tsv_mapping("mnxr_to_rhea.tsv")

    gene_count = annotate_genes(model)
    rxn_counts = annotate_reactions(model, seed_rhea=seed_rhea, mnxr_rhea=mnxr_rhea)

    summary = {
        "genes_ncbiprotein": gene_count,
        "reactions_brenda": rxn_counts["brenda"],
        "reactions_rhea": rxn_counts["rhea"],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(output_path))

    logger.info(
        f"{model_path.name}: +{gene_count} ncbiprotein genes, "
        f"+{rxn_counts['brenda']} brenda, +{rxn_counts['rhea']} rhea reactions"
    )
    return summary


def enrich_models_directory(
    models_dir: Path,
    output_dir: Path,
) -> Path:
    """Enrich all SBML models in a directory. Returns the output directory.

    Parameters
    ----------
    models_dir : Path
        Directory containing .xml model files (searched recursively).
    output_dir : Path
        Directory to write enriched models (flat structure).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    xml_files = list(models_dir.rglob("*.xml"))

    if not xml_files:
        logger.warning(f"No XML models found in {models_dir} — skipping annotation.")
        return output_dir

    logger.info(f"Enriching annotations for {len(xml_files)} model(s)...")

    # Load mappings once for all models
    seed_rhea = _load_tsv_mapping("seed_to_rhea.tsv")
    mnxr_rhea = _load_tsv_mapping("mnxr_to_rhea.tsv")

    totals = {"genes_ncbiprotein": 0, "reactions_brenda": 0, "reactions_rhea": 0}

    n = len(xml_files)
    for i, xml in enumerate(xml_files, 1):
        out_file = output_dir / xml.name
        model = cobra.io.read_sbml_model(str(xml))

        gene_count = annotate_genes(model)
        rxn_counts = annotate_reactions(model, seed_rhea=seed_rhea, mnxr_rhea=mnxr_rhea)

        cobra.io.write_sbml_model(model, str(out_file))

        totals["genes_ncbiprotein"] += gene_count
        totals["reactions_brenda"] += rxn_counts["brenda"]
        totals["reactions_rhea"] += rxn_counts["rhea"]

        logger.info(
            f"  [{i}/{n}] {xml.name}: +{gene_count} ncbiprotein, "
            f"+{rxn_counts['brenda']} brenda, +{rxn_counts['rhea']} rhea"
        )

    logger.info(
        f"Annotation enrichment complete: "
        f"+{totals['genes_ncbiprotein']} ncbiprotein genes, "
        f"+{totals['reactions_brenda']} brenda, "
        f"+{totals['reactions_rhea']} rhea reactions (across {len(xml_files)} models)"
    )

    return output_dir
