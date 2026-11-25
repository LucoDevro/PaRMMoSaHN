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
logging.getLogger("cobra.medium.boundary_types").setLevel(logging.WARNING)

# Note: bare IDs in cobra's `annotation` dict (e.g. "glc__D") are wrapped into
# MIRIAM identifiers.org URIs (e.g. https://identifiers.org/bigg.metabolite/glc__D)
# by cobra.io.write_sbml_model at write time. We deliberately store bare IDs in
# the dict — this is cobrapy's convention and produces the canonical URI form
# in the on-disk SBML. See tests/test_annotation.py for the regression test.

# SBO terms for SBML reaction/metabolite classification (MIRIAM standard)
_SBO = {
    "exchange": "SBO:0000627",  # exchange reaction
    "demand": "SBO:0000628",  # demand reaction
    "sink": "SBO:0000632",  # sink reaction
    "biomass": "SBO:0000629",  # biomass production
    "biochemical": "SBO:0000176",  # biochemical or transport reaction (default)
    "metabolite": "SBO:0000247",  # simple chemical
}

# Pattern for gapseq-mangled gene IDs:
#   {accession}_{version}_{start}_{end}  e.g. AEP93129_1_1_308
# We recover {accession}.{version} as the NCBI protein accession.
# Handles underscore-containing prefixes like NP_391935_1_1_308 → NP_391935.1
_GENE_ID_RE = re.compile(r"^(?P<acc>[A-Za-z]{2,3}_?\d+)_(?P<ver>\d+)_\d+_\d+$")

# Exchange reaction ID pattern for gapseq/SEED models: EX_cpd#####_e0
_SEED_EXCHANGE_RE = re.compile(r"^EX_(cpd\d+)_e0$")

# Metabolite ID pattern for gapseq/SEED models: cpd#####_{compartment}
# Matches cytoplasm (c0), external (e0), periplasm (p0), etc.
_SEED_MET_RE = re.compile(r"^(cpd\d+)_[a-z]\d*$")

# Internal reaction ID pattern: rxn#####_{compartment}
_SEED_RXN_RE = re.compile(r"^(rxn\d+)_[a-z]\d*$")

# SEED -> BiGG compartment-suffix translation. SEED uses {letter}{digit}
# (e.g. c0, e0, p0); BiGG drops the digit. We expose this as a separate
# constant so rename_to_bigg's behaviour is greppable.
_SEED_TO_BIGG_COMPARTMENT = {"c0": "c", "e0": "e", "p0": "p"}

# Pattern for Bakta-style gene IDs produced by PPanGGOLiN + gapseq:
#   G_gnl_Bakta_{assembly}_{ver}_{locus}_{start}_{end}
#   e.g. G_gnl_Bakta_GCF_000144625_1_00215_3_298
_BAKTA_GENE_ID_RE = re.compile(
    r"^(?:G_)?gnl_Bakta_(?P<core>.+)_(?P<start>\d+)_(?P<end>\d+)$"
)


def _core_to_bakta_protein_id(core: str) -> str:
    """Reconstruct a GenBank protein_id from the underscore-flattened core.

    ``GCF_000144625_1_00215`` → ``gnl|Bakta|GCF_000144625.1_00215``

    The assembly accession follows the pattern ``GCF_NNNNNNNNN_V`` where
    V is the version digit.  We restore the ``.`` between the accession
    number and the version, then reattach the locus suffix.
    """
    # Split: ['GCF', '000144625', '1', '00215']
    parts = core.split("_")
    if len(parts) >= 4:
        # GCF_000144625.1_00215  (prefix_number.version_locus)
        accession = f"{parts[0]}_{parts[1]}.{parts[2]}"
        locus = "_".join(parts[3:])
        return f"gnl|Bakta|{accession}_{locus}"
    return f"gnl|Bakta|{core}"


def load_gene_mapping(tsv_path: Path) -> dict[str, str]:
    """Load a protein_id → RefSeq accession mapping from a two-column TSV."""
    mapping: dict[str, str] = {}
    if not tsv_path.exists():
        return mapping
    for line in tsv_path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t", 1)
        if len(parts) == 2:
            mapping[parts[0]] = parts[1]
    return mapping


def _load_tsv_mapping(tsv_name: str) -> dict[str, str]:
    """Load a two-column TSV from the bundled data directory into a dict.

    Returns an empty mapping (with a warning) if the bundled resource is missing
    or empty, so a packaging gap degrades enrichment gracefully rather than
    aborting the whole run with ``FileNotFoundError`` / ``StopIteration``.
    """
    resource = files("parmmosahn.data").joinpath(tsv_name)
    mapping: dict[str, str] = {}
    if not resource.is_file():
        logger.warning(
            f"Bundled mapping {tsv_name} not found — enrichment from this source "
            f"will be skipped. Run `parmmosahn doctor` to verify the install."
        )
        return mapping
    with resource.open("r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # skip header, tolerate an empty file
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


def annotate_genes(
    model: cobra.Model, gene_mapping: dict[str, str] | None = None
) -> int:
    """Recover NCBI protein accessions from gapseq-mangled gene IDs.

    Supports two ID formats:
    - NCBI-style: ``AEP93129_1_1_308`` → ``AEP93129.1``
    - Bakta-style: ``G_gnl_Bakta_GCF_000144625_1_00215_3_298``
      → looked up via *gene_mapping* dict (protein_id → RefSeq accession)

    Returns the number of genes annotated.
    """
    count = 0
    for gene in model.genes:
        if _get_annotation(gene, "ncbiprotein"):
            continue
        # Try NCBI-style regex first
        m = _GENE_ID_RE.match(gene.id)
        if m:
            accession = f"{m.group('acc')}.{m.group('ver')}"
            if _add_annotation(gene, "ncbiprotein", accession):
                count += 1
            continue
        # Try Bakta-style with mapping lookup
        if gene_mapping:
            bm = _BAKTA_GENE_ID_RE.match(gene.id)
            if bm:
                bakta_id = _core_to_bakta_protein_id(bm.group("core"))
                ncbi_acc = gene_mapping.get(bakta_id)
                if ncbi_acc and _add_annotation(gene, "ncbiprotein", ncbi_acc):
                    count += 1

    if count == 0 and any("Bakta" in g.id or "gnl" in g.id for g in model.genes):
        logger.warning(
            "No NCBI protein accessions recovered. Gene IDs appear to be "
            "Bakta-style. Ensure .gene_mapping.tsv files exist in the "
            "proteomes/ directory (generated from GenBank db_xref fields)."
        )
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


def annotate_exchanges(
    model: cobra.Model,
    seed_bigg: dict[str, str] | None = None,
) -> int:
    """Add BiGG annotations to exchange reactions and their boundary metabolites.

    For each SEED-format exchange reaction (``EX_cpd#####_e0``):
    - ``bigg.metabolite`` is written on the boundary metabolite
    - ``bigg.reaction``  is written on the reaction itself (e.g. ``EX_glc__D_e``)

    Uses the bundled ``seed_to_bigg.tsv`` (2900+ entries built from MetaNetX
    chem_xref plus supervisor-verified overrides) unless *seed_bigg* is supplied.

    Returns the number of exchange reactions annotated with ``bigg.reaction``
    (idempotent: already-present annotations are not re-added and not counted).
    """
    if seed_bigg is None:
        seed_bigg = _load_tsv_mapping("seed_to_bigg.tsv")

    count = 0
    for rxn in model.exchanges:
        m = _SEED_EXCHANGE_RE.match(rxn.id)
        if not m:
            continue
        bigg_id = seed_bigg.get(m.group(1))
        if not bigg_id:
            continue
        for met in rxn.metabolites:
            _add_annotation(met, "bigg.metabolite", bigg_id)
        if _add_annotation(rxn, "bigg.reaction", f"EX_{bigg_id}_e"):
            count += 1

    return count


def annotate_metabolites(
    model: cobra.Model,
    seed_bigg: dict[str, str] | None = None,
) -> int:
    """Add ``bigg.metabolite`` to every SEED-format metabolite in the model.

    Iterates ``model.metabolites`` (all compartments: cytoplasm, external,
    periplasm, …).  Boundary metabolites of exchange reactions already
    annotated by :func:`annotate_exchanges` are skipped (idempotent).

    Returns the number of metabolites newly annotated.
    """
    if seed_bigg is None:
        seed_bigg = _load_tsv_mapping("seed_to_bigg.tsv")

    count = 0
    for met in model.metabolites:
        m = _SEED_MET_RE.match(met.id)
        if not m:
            continue
        bigg_id = seed_bigg.get(m.group(1))
        if not bigg_id:
            continue
        if _add_annotation(met, "bigg.metabolite", bigg_id):
            count += 1

    return count


def annotate_sbo_terms(model: cobra.Model) -> dict[str, int]:
    """Add SBO terms to every reaction and metabolite that lacks one.

    Classification for reactions:
    - Exchange reactions (``model.exchanges``)  → SBO:0000627
    - Demand reactions   (``model.demands``)    → SBO:0000628
    - Sink reactions     (``model.sinks``)      → SBO:0000632
    - Reactions with "biomass" in their ID/name → SBO:0000629
    - All others                                → SBO:0000176

    All metabolites receive SBO:0000247 ("simple chemical").

    Returns ``{"reactions": n_added, "metabolites": n_added}``.
    """
    exchange_ids = {rxn.id for rxn in model.exchanges}
    demand_ids = {rxn.id for rxn in model.demands}
    sink_ids = {rxn.id for rxn in model.sinks}

    rxn_count = 0
    for rxn in model.reactions:
        if rxn.id in exchange_ids:
            sbo = _SBO["exchange"]
        elif rxn.id in demand_ids:
            sbo = _SBO["demand"]
        elif rxn.id in sink_ids:
            sbo = _SBO["sink"]
        elif "biomass" in rxn.id.lower() or "biomass" in (rxn.name or "").lower():
            sbo = _SBO["biomass"]
        else:
            sbo = _SBO["biochemical"]
        if _add_annotation(rxn, "sbo", sbo):
            rxn_count += 1

    met_count = sum(
        1
        for met in model.metabolites
        if _add_annotation(met, "sbo", _SBO["metabolite"])
    )

    return {"reactions": rxn_count, "metabolites": met_count}


def _pick_bigg_reaction_id(
    rxn: "cobra.Reaction", seed_rxn_bigg: dict[str, str]
) -> str | None:
    """Policy B: prefer an existing bigg.reaction annotation over the TSV.

    Gapseq sometimes emits a list of compartment-pinned BiGG variants on a
    single SEED reaction (e.g. ['PPA', 'PPA_1pp', 'IDPh', 'IDPm']). When such
    a list is present we take the first entry — it has compartment context
    that the MNX chain in the bundled TSV lacks. Only when no bigg.reaction
    annotation exists at all do we fall back to seed_rxn_bigg.
    """
    existing = (rxn.annotation or {}).get("bigg.reaction")
    if existing:
        if isinstance(existing, list):
            return existing[0] if existing else None
        return existing
    m = _SEED_RXN_RE.match(rxn.id)
    if not m:
        return None
    return seed_rxn_bigg.get(m.group(1))


def rename_to_bigg(
    model: cobra.Model,
    seed_bigg: dict[str, str] | None = None,
    seed_rxn_bigg: dict[str, str] | None = None,
) -> dict[str, int]:
    """Rewrite SEED-namespace primary IDs to BiGG IDs in-place (opt-in).

    Three rename classes:
    - Metabolites ``cpd#####_{c0,e0,p0}`` -> ``{bigg}_{c,e,p}`` via
      ``seed_to_bigg.tsv`` (~52 % coverage on the audit cohort).
    - Exchange reactions ``EX_cpd#####_e0`` -> ``EX_{bigg}_e``, driven by the
      same metabolite mapping.
    - Internal reactions ``rxn#####_{c0,...}`` -> BiGG ID. Policy B (see
      :func:`_pick_bigg_reaction_id`): prefer any pre-existing ``bigg.reaction``
      annotation (gapseq writes these with compartment context); fall back to
      ``seed_rxn_to_bigg.tsv``.

    The original SEED ID is already carried in ``seed.compound`` /
    ``seed.reaction`` for >99 % / >91 % of entities post-gapseq, so renaming
    does not lose provenance. Collisions with an ID already in the model are
    skipped with a WARN log; running on an already-renamed model is a no-op
    (the regex match fails on BiGG-form IDs).

    Returns ``{'metabolites', 'reactions_internal', 'reactions_exchange',
    'collisions_skipped'}``.
    """
    if seed_bigg is None:
        seed_bigg = _load_tsv_mapping("seed_to_bigg.tsv")
    if seed_rxn_bigg is None:
        seed_rxn_bigg = _load_tsv_mapping("seed_rxn_to_bigg.tsv")

    counts = {
        "metabolites": 0,
        "reactions_internal": 0,
        "reactions_exchange": 0,
        "collisions_skipped": 0,
    }

    taken_met_ids = {met.id for met in model.metabolites}
    for met in list(model.metabolites):
        mt = _SEED_MET_RE.match(met.id)
        if not mt:
            continue
        cpd, seed_comp = mt.group(1), met.id[len(mt.group(1)) + 1 :]
        bigg = seed_bigg.get(cpd)
        if not bigg:
            continue
        new_comp = _SEED_TO_BIGG_COMPARTMENT.get(seed_comp, seed_comp)
        new_id = f"{bigg}_{new_comp}"
        if new_id == met.id:
            continue
        if new_id in taken_met_ids:
            counts["collisions_skipped"] += 1
            logger.warning(
                "rename_to_bigg: metabolite %s -> %s collides; keeping SEED ID",
                met.id,
                new_id,
            )
            continue
        taken_met_ids.discard(met.id)
        taken_met_ids.add(new_id)
        met.id = new_id
        counts["metabolites"] += 1

    exchange_id_set = {ex.id for ex in model.exchanges}
    taken_rxn_ids = {rxn.id for rxn in model.reactions}
    for rxn in list(model.reactions):
        if rxn.id not in exchange_id_set:
            continue
        em = _SEED_EXCHANGE_RE.match(rxn.id)
        if not em:
            continue
        bigg = seed_bigg.get(em.group(1))
        if not bigg:
            continue
        new_id = f"EX_{bigg}_e"
        if new_id == rxn.id:
            continue
        if new_id in taken_rxn_ids:
            counts["collisions_skipped"] += 1
            logger.warning(
                "rename_to_bigg: exchange %s -> %s collides; keeping SEED ID",
                rxn.id,
                new_id,
            )
            continue
        taken_rxn_ids.discard(rxn.id)
        taken_rxn_ids.add(new_id)
        rxn.id = new_id
        counts["reactions_exchange"] += 1

    for rxn in list(model.reactions):
        # _SEED_RXN_RE only matches `rxn\d+_x\d*`, so non-SEED reactions
        # (biomass, demand, sink, exchanges, already-renamed) are skipped.
        if not _SEED_RXN_RE.match(rxn.id):
            continue
        picked = _pick_bigg_reaction_id(rxn, seed_rxn_bigg)
        if not picked or picked == rxn.id:
            continue
        if picked in taken_rxn_ids:
            counts["collisions_skipped"] += 1
            logger.warning(
                "rename_to_bigg: reaction %s -> %s collides; keeping SEED ID",
                rxn.id,
                picked,
            )
            continue
        taken_rxn_ids.discard(rxn.id)
        taken_rxn_ids.add(picked)
        rxn.id = picked
        counts["reactions_internal"] += 1

    return counts


def polish_models(
    models_dir: Path,
    output_dir: Path,
    proteomes_dir: Path | None = None,
    rename: bool = False,
    skip_existing: bool = True,
) -> Path:
    """Enrich all SBML models in a directory. Returns the output directory.

    With *skip_existing* (default), a model whose enriched output already exists
    and is non-empty is skipped, so a resumed run does not re-enrich it.

    Parameters
    ----------
    models_dir : Path
        Directory containing .xml model files (searched recursively).
    output_dir : Path
        Directory to write enriched models (flat structure).
    proteomes_dir : Path, optional
        Directory containing ``*.gene_mapping.tsv`` files (one per strain).
        If provided, Bakta-style gene IDs are resolved to NCBI accessions.
    rename : bool
        If True, also rewrite SEED-namespace primary IDs to BiGG IDs via
        :func:`rename_to_bigg` after annotation enrichment. Opt-in.
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
    seed_bigg = _load_tsv_mapping("seed_to_bigg.tsv")
    seed_rxn_bigg = _load_tsv_mapping("seed_rxn_to_bigg.tsv") if rename else {}

    totals: dict[str, int] = {
        "genes_ncbiprotein": 0,
        "reactions_brenda": 0,
        "reactions_rhea": 0,
        "exchanges_bigg": 0,
        "metabolites_bigg": 0,
        "sbo_reactions": 0,
        "sbo_metabolites": 0,
        "renamed_metabolites": 0,
        "renamed_reactions_internal": 0,
        "renamed_reactions_exchange": 0,
        "rename_collisions_skipped": 0,
    }

    n = len(xml_files)
    failed: list[str] = []
    skipped = 0
    for i, xml in enumerate(xml_files, 1):
        out_file = output_dir / xml.name
        if skip_existing and out_file.exists() and out_file.stat().st_size > 0:
            skipped += 1
            continue
        try:
            model = cobra.io.read_sbml_model(str(xml))

            # Load per-strain gene mapping if available
            gene_mapping: dict[str, str] | None = None
            if proteomes_dir:
                mapping_tsv = proteomes_dir / f"{xml.stem}.gene_mapping.tsv"
                if mapping_tsv.exists():
                    gene_mapping = load_gene_mapping(mapping_tsv)

            gene_count = annotate_genes(model, gene_mapping=gene_mapping)
            rxn_counts = annotate_reactions(
                model, seed_rhea=seed_rhea, mnxr_rhea=mnxr_rhea
            )
            exch_count = annotate_exchanges(model, seed_bigg=seed_bigg)
            met_count = annotate_metabolites(model, seed_bigg=seed_bigg)
            sbo_counts = annotate_sbo_terms(model)
            if rename:
                rename_counts = rename_to_bigg(
                    model, seed_bigg=seed_bigg, seed_rxn_bigg=seed_rxn_bigg
                )
            else:
                rename_counts = {
                    "metabolites": 0,
                    "reactions_internal": 0,
                    "reactions_exchange": 0,
                    "collisions_skipped": 0,
                }

            cobra.io.write_sbml_model(model, str(out_file))
        except Exception as exc:  # one bad model must not abort the batch
            failed.append(xml.name)
            logger.error(f"Annotation enrichment failed for {xml.name}: {exc}")
            continue

        totals["genes_ncbiprotein"] += gene_count
        totals["reactions_brenda"] += rxn_counts["brenda"]
        totals["reactions_rhea"] += rxn_counts["rhea"]
        totals["exchanges_bigg"] += exch_count
        totals["metabolites_bigg"] += met_count
        totals["sbo_reactions"] += sbo_counts["reactions"]
        totals["sbo_metabolites"] += sbo_counts["metabolites"]
        totals["renamed_metabolites"] += rename_counts["metabolites"]
        totals["renamed_reactions_internal"] += rename_counts["reactions_internal"]
        totals["renamed_reactions_exchange"] += rename_counts["reactions_exchange"]
        totals["rename_collisions_skipped"] += rename_counts["collisions_skipped"]

        logger.debug(
            f"  [{i}/{n}] {xml.name}: +{gene_count} ncbiprotein, "
            f"+{rxn_counts['brenda']} brenda, +{rxn_counts['rhea']} rhea, "
            f"+{exch_count} bigg.reaction, +{met_count} bigg.metabolite, "
            f"+{sbo_counts['reactions']} sbo rxns, +{sbo_counts['metabolites']} sbo mets"
        )

    if skipped:
        logger.info(
            f"Annotation enrichment skipped {skipped}/{n} model(s) already enriched."
        )
    if failed:
        logger.warning(
            f"Annotation enrichment skipped {len(failed)}/{n} model(s) due to "
            f"errors: {', '.join(failed)}"
        )

    logger.info(
        f"Annotation enrichment complete: "
        f"+{totals['genes_ncbiprotein']} ncbiprotein genes, "
        f"+{totals['reactions_brenda']} brenda, "
        f"+{totals['reactions_rhea']} rhea, "
        f"+{totals['exchanges_bigg']} bigg.reaction exchanges, "
        f"+{totals['metabolites_bigg']} bigg.metabolite internal, "
        f"+{totals['sbo_reactions']} sbo reactions, "
        f"+{totals['sbo_metabolites']} sbo metabolites "
        f"(across {len(xml_files)} models)"
    )

    return output_dir
