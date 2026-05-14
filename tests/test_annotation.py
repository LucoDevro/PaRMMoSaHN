from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from parmmosahn.annotation import (
    _add_annotation,
    _get_annotation,
    annotate_exchanges,
    annotate_genes,
    annotate_metabolites,
    annotate_reactions,
    annotate_sbo_terms,
    enrich_model_annotations,
    rename_to_bigg,
)
from parmmosahn.cli import cli


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_gene(gene_id: str, annotation: dict | None = None) -> MagicMock:
    g = MagicMock()
    g.id = gene_id
    g.annotation = annotation if annotation is not None else {}
    return g


def _make_reaction(rxn_id: str, annotation: dict | None = None) -> MagicMock:
    r = MagicMock()
    r.id = rxn_id
    r.annotation = annotation if annotation is not None else {}
    return r


def _make_model(genes=None, reactions=None) -> MagicMock:
    model = MagicMock()
    model.genes = genes or []
    model.reactions = reactions or []
    return model


# ---------------------------------------------------------------------------
# Unit tests: _get_annotation / _add_annotation
# ---------------------------------------------------------------------------
def test_get_annotation_missing_key():
    entity = MagicMock()
    entity.annotation = {"ec-code": "1.1.1.1"}
    assert _get_annotation(entity, "brenda") == []


def test_get_annotation_string_value():
    entity = MagicMock()
    entity.annotation = {"ec-code": "1.1.1.1"}
    assert _get_annotation(entity, "ec-code") == ["1.1.1.1"]


def test_get_annotation_list_value():
    entity = MagicMock()
    entity.annotation = {"ec-code": ["1.1.1.1", "2.2.2.2"]}
    assert _get_annotation(entity, "ec-code") == ["1.1.1.1", "2.2.2.2"]


def test_add_annotation_new_key():
    entity = MagicMock()
    entity.annotation = {}
    assert _add_annotation(entity, "brenda", "1.1.1.1") is True
    assert entity.annotation["brenda"] == "1.1.1.1"


def test_add_annotation_duplicate_string():
    entity = MagicMock()
    entity.annotation = {"brenda": "1.1.1.1"}
    assert _add_annotation(entity, "brenda", "1.1.1.1") is False


def test_add_annotation_promotes_to_list():
    entity = MagicMock()
    entity.annotation = {"brenda": "1.1.1.1"}
    assert _add_annotation(entity, "brenda", "2.2.2.2") is True
    assert entity.annotation["brenda"] == ["1.1.1.1", "2.2.2.2"]


def test_add_annotation_appends_to_list():
    entity = MagicMock()
    entity.annotation = {"brenda": ["1.1.1.1"]}
    assert _add_annotation(entity, "brenda", "2.2.2.2") is True
    assert "2.2.2.2" in entity.annotation["brenda"]


def test_add_annotation_duplicate_in_list():
    entity = MagicMock()
    entity.annotation = {"brenda": ["1.1.1.1"]}
    assert _add_annotation(entity, "brenda", "1.1.1.1") is False


# ---------------------------------------------------------------------------
# Unit tests: annotate_genes
# ---------------------------------------------------------------------------
def test_annotate_genes_recovers_accessions():
    """Gene IDs like AEP93129_1_1_308 → ncbiprotein: AEP93129.1"""
    genes = [
        _make_gene("AEP93129_1_1_308"),
        _make_gene("NP_391935_1_1_308"),
        _make_gene("WP_013308093_1_1_308"),
    ]
    model = _make_model(genes=genes)

    count = annotate_genes(model)

    assert count == 3
    assert genes[0].annotation["ncbiprotein"] == "AEP93129.1"
    assert genes[1].annotation["ncbiprotein"] == "NP_391935.1"
    assert genes[2].annotation["ncbiprotein"] == "WP_013308093.1"


def test_annotate_genes_skips_existing():
    """Genes already annotated with ncbiprotein should be left alone."""
    genes = [
        _make_gene("AEP93129_1_1_308", {"ncbiprotein": "AEP93129.1"}),
    ]
    model = _make_model(genes=genes)

    count = annotate_genes(model)
    assert count == 0


def test_annotate_genes_bakta_style_with_G_prefix():
    """Bakta-style gene IDs with G_ prefix should be resolved via gene_mapping."""
    genes = [
        _make_gene("G_gnl_Bakta_GCF_000144625_1_00215_3_298"),
    ]
    model = _make_model(genes=genes)
    gene_mapping = {"gnl|Bakta|GCF_000144625.1_00215": "WP_012345678.1"}

    count = annotate_genes(model, gene_mapping=gene_mapping)

    assert count == 1
    assert genes[0].annotation["ncbiprotein"] == "WP_012345678.1"


def test_annotate_genes_bakta_style_without_G_prefix():
    """Bakta-style gene IDs without G_ prefix (real Clostridiales output)."""
    genes = [
        _make_gene("gnl_Bakta_GCF_000007625_1_00510_3_298"),
        _make_gene("gnl_Bakta_GCF_000007625_1_05095_3_592"),
    ]
    model = _make_model(genes=genes)
    gene_mapping = {
        "gnl|Bakta|GCF_000007625.1_00510": "WP_035110684.1",
        "gnl|Bakta|GCF_000007625.1_05095": "WP_011098500.1",
    }

    count = annotate_genes(model, gene_mapping=gene_mapping)

    assert count == 2
    assert genes[0].annotation["ncbiprotein"] == "WP_035110684.1"
    assert genes[1].annotation["ncbiprotein"] == "WP_011098500.1"


def test_annotate_genes_bakta_no_mapping_warns(caplog):
    """Bakta-style IDs without gene_mapping should warn (no crash)."""
    genes = [_make_gene("gnl_Bakta_GCF_000007625_1_00510_3_298")]
    model = _make_model(genes=genes)

    count = annotate_genes(model, gene_mapping=None)

    assert count == 0
    assert "Bakta-style" in caplog.text


def test_annotate_genes_unmatched_id():
    """Gene IDs that don't match the pattern should be skipped."""
    genes = [
        _make_gene("some_random_id"),
        _make_gene("s0001"),
    ]
    model = _make_model(genes=genes)

    count = annotate_genes(model)
    assert count == 0


# ---------------------------------------------------------------------------
# Unit tests: annotate_reactions
# ---------------------------------------------------------------------------
def test_annotate_reactions_brenda_from_ec():
    """BRENDA annotation should be added from ec-code."""
    rxns = [
        _make_reaction("rxn1", {"ec-code": "3.6.1.1"}),
        _make_reaction("rxn2", {"ec-code": ["1.1.1.1", "2.2.2.2"]}),
        _make_reaction("rxn3", {}),  # no ec-code
    ]
    model = _make_model(reactions=rxns)

    counts = annotate_reactions(model, seed_rhea={}, mnxr_rhea={})

    assert counts["brenda"] == 2
    assert rxns[0].annotation["brenda"] == "3.6.1.1"
    assert rxns[1].annotation["brenda"] == "1.1.1.1"


def test_annotate_reactions_brenda_skips_existing():
    """Don't overwrite existing BRENDA annotations."""
    rxns = [
        _make_reaction("rxn1", {"ec-code": "3.6.1.1", "brenda": "3.6.1.1"}),
    ]
    model = _make_model(reactions=rxns)

    counts = annotate_reactions(model, seed_rhea={}, mnxr_rhea={})
    assert counts["brenda"] == 0


def test_annotate_reactions_rhea_from_seed():
    """RHEA annotation via seed.reaction → bundled mapping."""
    seed_rhea = {"rxn00001": "24576"}
    rxns = [
        _make_reaction("rxn1", {"seed.reaction": "rxn00001"}),
    ]
    model = _make_model(reactions=rxns)

    counts = annotate_reactions(model, seed_rhea=seed_rhea, mnxr_rhea={})

    assert counts["rhea"] == 1
    assert rxns[0].annotation["rhea"] == "24576"


def test_annotate_reactions_rhea_from_metanetx():
    """RHEA annotation via metanetx.reaction when seed.reaction has no mapping."""
    mnxr_rhea = {"MNXR03": "34979"}
    rxns = [
        _make_reaction("rxn1", {"metanetx.reaction": "MNXR03"}),
    ]
    model = _make_model(reactions=rxns)

    counts = annotate_reactions(model, seed_rhea={}, mnxr_rhea=mnxr_rhea)

    assert counts["rhea"] == 1
    assert rxns[0].annotation["rhea"] == "34979"


def test_annotate_reactions_rhea_seed_preferred_over_metanetx():
    """seed.reaction mapping should be tried first."""
    seed_rhea = {"rxn00001": "24576"}
    mnxr_rhea = {"MNXR03": "34979"}
    rxns = [
        _make_reaction(
            "rxn1",
            {"seed.reaction": "rxn00001", "metanetx.reaction": "MNXR03"},
        ),
    ]
    model = _make_model(reactions=rxns)

    counts = annotate_reactions(model, seed_rhea=seed_rhea, mnxr_rhea=mnxr_rhea)

    assert counts["rhea"] == 1
    assert rxns[0].annotation["rhea"] == "24576"  # from seed, not metanetx


def test_annotate_reactions_rhea_skips_existing():
    """Don't overwrite existing RHEA annotations."""
    seed_rhea = {"rxn00001": "24576"}
    rxns = [
        _make_reaction("rxn1", {"seed.reaction": "rxn00001", "rhea": "99999"}),
    ]
    model = _make_model(reactions=rxns)

    counts = annotate_reactions(model, seed_rhea=seed_rhea, mnxr_rhea={})
    assert counts["rhea"] == 0


# ---------------------------------------------------------------------------
# Unit tests: annotate_exchanges
# ---------------------------------------------------------------------------
def _make_metabolite(met_id: str, annotation: dict | None = None) -> MagicMock:
    met = MagicMock()
    met.id = met_id
    met.annotation = annotation if annotation is not None else {}
    return met


def _make_exchange_reaction(
    rxn_id: str, met_id: str, met_annotation: dict | None = None
):
    """Return (reaction_mock, metabolite_mock) for an exchange reaction."""
    met = MagicMock()
    met.id = met_id
    met.annotation = met_annotation if met_annotation is not None else {}
    rxn = MagicMock()
    rxn.id = rxn_id
    rxn.annotation = {}
    rxn.metabolites = {met: -1.0}
    return rxn, met


def test_annotate_exchanges_adds_bigg_annotation():
    """EX_cpd00027_e0 → bigg.metabolite on boundary met + bigg.reaction on reaction."""
    seed_bigg = {"cpd00027": "glc__D", "cpd00082": "fru"}
    rxn1, met1 = _make_exchange_reaction("EX_cpd00027_e0", "cpd00027_e0")
    rxn2, met2 = _make_exchange_reaction("EX_cpd00082_e0", "cpd00082_e0")
    model = MagicMock()
    model.exchanges = [rxn1, rxn2]

    count = annotate_exchanges(model, seed_bigg=seed_bigg)

    assert count == 2
    assert met1.annotation["bigg.metabolite"] == "glc__D"
    assert rxn1.annotation["bigg.reaction"] == "EX_glc__D_e"
    assert met2.annotation["bigg.metabolite"] == "fru"
    assert rxn2.annotation["bigg.reaction"] == "EX_fru_e"


def test_annotate_exchanges_skips_unmapped_seed_cpd():
    """SEED IDs not in seed_bigg are silently skipped."""
    seed_bigg = {"cpd00027": "glc__D"}
    rxn, met = _make_exchange_reaction("EX_cpd99999_e0", "cpd99999_e0")
    model = MagicMock()
    model.exchanges = [rxn]

    count = annotate_exchanges(model, seed_bigg=seed_bigg)

    assert count == 0
    assert "bigg.metabolite" not in met.annotation
    assert "bigg.reaction" not in rxn.annotation


def test_annotate_exchanges_skips_non_seed_exchange_ids():
    """BiGG-style exchange IDs (from ModelPolisher) don't match the SEED pattern."""
    seed_bigg = {"cpd00027": "glc__D"}
    rxn, met = _make_exchange_reaction("EX_glc__D_e", "glc__D_e")
    model = MagicMock()
    model.exchanges = [rxn]

    count = annotate_exchanges(model, seed_bigg=seed_bigg)

    assert count == 0
    assert "bigg.reaction" not in rxn.annotation


def test_annotate_exchanges_skips_existing_bigg_reaction():
    """Exchange reactions already annotated with bigg.reaction are not re-added."""
    seed_bigg = {"cpd00027": "glc__D"}
    rxn, met = _make_exchange_reaction("EX_cpd00027_e0", "cpd00027_e0")
    rxn.annotation = {"bigg.reaction": "EX_glc__D_e"}
    model = MagicMock()
    model.exchanges = [rxn]

    count = annotate_exchanges(model, seed_bigg=seed_bigg)

    assert count == 0  # bigg.reaction already present


def test_annotate_exchanges_empty_exchange_list():
    """A model with no exchange reactions returns 0."""
    model = MagicMock()
    model.exchanges = []

    count = annotate_exchanges(model, seed_bigg={"cpd00027": "glc__D"})
    assert count == 0


# ---------------------------------------------------------------------------
# Unit tests: annotate_metabolites
# ---------------------------------------------------------------------------
def test_annotate_metabolites_adds_bigg_internal():
    """Intracellular SEED metabolites (cpd#####_c0) get bigg.metabolite."""
    seed_bigg = {"cpd00027": "glc__D", "cpd00004": "nadh"}
    mets = [
        _make_metabolite("cpd00027_c0"),
        _make_metabolite("cpd00004_c0"),
    ]
    model = MagicMock()
    model.metabolites = mets

    count = annotate_metabolites(model, seed_bigg=seed_bigg)

    assert count == 2
    assert mets[0].annotation["bigg.metabolite"] == "glc__D"
    assert mets[1].annotation["bigg.metabolite"] == "nadh"


def test_annotate_metabolites_all_compartments():
    """Matches cytoplasm (_c0), external (_e0), and periplasm (_p0)."""
    seed_bigg = {"cpd00027": "glc__D"}
    mets = [
        _make_metabolite("cpd00027_c0"),
        _make_metabolite("cpd00027_e0"),
        _make_metabolite("cpd00027_p0"),
    ]
    model = MagicMock()
    model.metabolites = mets

    count = annotate_metabolites(model, seed_bigg=seed_bigg)
    assert count == 3


def test_annotate_metabolites_skips_non_seed_ids():
    """BiGG-style metabolite IDs (from ModelPolisher) are not matched."""
    seed_bigg = {"cpd00027": "glc__D"}
    mets = [_make_metabolite("glc__D_c"), _make_metabolite("glc__D_e")]
    model = MagicMock()
    model.metabolites = mets

    count = annotate_metabolites(model, seed_bigg=seed_bigg)
    assert count == 0


def test_annotate_metabolites_skips_unmapped():
    """SEED metabolites not in seed_bigg are silently skipped."""
    seed_bigg = {"cpd00027": "glc__D"}
    mets = [_make_metabolite("cpd99999_c0")]
    model = MagicMock()
    model.metabolites = mets

    count = annotate_metabolites(model, seed_bigg=seed_bigg)
    assert count == 0


def test_annotate_metabolites_skips_already_annotated():
    """Metabolites already carrying bigg.metabolite are not re-annotated."""
    seed_bigg = {"cpd00027": "glc__D"}
    mets = [_make_metabolite("cpd00027_c0", {"bigg.metabolite": "glc__D"})]
    model = MagicMock()
    model.metabolites = mets

    count = annotate_metabolites(model, seed_bigg=seed_bigg)
    assert count == 0


def test_annotate_metabolites_idempotent_after_annotate_exchanges():
    """Boundary metabolites annotated by annotate_exchanges are skipped by annotate_metabolites."""
    seed_bigg = {"cpd00027": "glc__D"}
    rxn, boundary_met = _make_exchange_reaction("EX_cpd00027_e0", "cpd00027_e0")
    model_ex = MagicMock()
    model_ex.exchanges = [rxn]
    # First pass: annotate_exchanges annotates boundary_met and rxn
    exch_count = annotate_exchanges(model_ex, seed_bigg=seed_bigg)
    assert exch_count == 1
    assert boundary_met.annotation["bigg.metabolite"] == "glc__D"

    # Second pass: annotate_metabolites skips boundary_met (already annotated)
    model_ex.metabolites = [boundary_met]
    met_count = annotate_metabolites(model_ex, seed_bigg=seed_bigg)
    assert met_count == 0


# ---------------------------------------------------------------------------
# Unit tests: annotate_sbo_terms
# ---------------------------------------------------------------------------
def _make_sbo_model(
    exchange_ids=(), demand_ids=(), sink_ids=(), extra_rxn_ids=(), met_ids=()
):
    """Minimal model mock for SBO tests."""

    def _rxn(rid, name=""):
        r = MagicMock()
        r.id = rid
        r.name = name
        r.annotation = {}
        return r

    exchanges = [_rxn(i) for i in exchange_ids]
    demands = [_rxn(i) for i in demand_ids]
    sinks = [_rxn(i) for i in sink_ids]
    extras = [_rxn(i) for i in extra_rxn_ids]
    mets = [_make_metabolite(m) for m in met_ids]

    model = MagicMock()
    model.exchanges = exchanges
    model.demands = demands
    model.sinks = sinks
    model.reactions = exchanges + demands + sinks + extras
    model.metabolites = mets
    return model, exchanges, demands, sinks, extras, mets


def test_annotate_sbo_exchange_reactions():
    """Exchange reactions get SBO:0000627."""
    model, exchanges, *_ = _make_sbo_model(
        exchange_ids=["EX_cpd00027_e0", "EX_cpd00082_e0"]
    )

    counts = annotate_sbo_terms(model)

    assert counts["reactions"] == 2
    assert exchanges[0].annotation["sbo"] == "SBO:0000627"
    assert exchanges[1].annotation["sbo"] == "SBO:0000627"


def test_annotate_sbo_demand_and_sink():
    """Demand → SBO:0000628, sink → SBO:0000632."""
    model, _, demands, sinks, *_ = _make_sbo_model(
        demand_ids=["DM_cpd00001_c0"], sink_ids=["SK_cpd00002_c0"]
    )

    annotate_sbo_terms(model)

    assert demands[0].annotation["sbo"] == "SBO:0000628"
    assert sinks[0].annotation["sbo"] == "SBO:0000632"


def test_annotate_sbo_biomass_by_id():
    """Reactions with 'biomass' in ID get SBO:0000629."""
    model, _, _, _, extras, _ = _make_sbo_model(extra_rxn_ids=["biomass_reaction"])

    annotate_sbo_terms(model)

    assert extras[0].annotation["sbo"] == "SBO:0000629"


def test_annotate_sbo_biomass_by_name():
    """Reactions with 'biomass' in name (not ID) also get SBO:0000629."""

    def _rxn_named(rid, name):
        r = MagicMock()
        r.id = rid
        r.name = name
        r.annotation = {}
        return r

    rxn = _rxn_named("OBJSYN", "Biomass production")
    model = MagicMock()
    model.exchanges = []
    model.demands = []
    model.sinks = []
    model.reactions = [rxn]
    model.metabolites = []

    annotate_sbo_terms(model)

    assert rxn.annotation["sbo"] == "SBO:0000629"


def test_annotate_sbo_biochemical_default():
    """Non-classified reactions fall back to SBO:0000176."""
    model, _, _, _, extras, _ = _make_sbo_model(extra_rxn_ids=["PGI", "PFK"])

    annotate_sbo_terms(model)

    assert extras[0].annotation["sbo"] == "SBO:0000176"
    assert extras[1].annotation["sbo"] == "SBO:0000176"


def test_annotate_sbo_metabolites():
    """All metabolites receive SBO:0000247."""
    model, *_, mets = _make_sbo_model(met_ids=["cpd00027_c0", "cpd00027_e0"])

    counts = annotate_sbo_terms(model)

    assert counts["metabolites"] == 2
    assert mets[0].annotation["sbo"] == "SBO:0000247"
    assert mets[1].annotation["sbo"] == "SBO:0000247"


def test_annotate_sbo_skips_existing():
    """Entities already carrying an sbo annotation are not overwritten."""
    model, _, _, _, extras, mets = _make_sbo_model(
        extra_rxn_ids=["PGI"], met_ids=["cpd00027_c0"]
    )
    extras[0].annotation["sbo"] = "SBO:0000176"
    mets[0].annotation["sbo"] = "SBO:0000247"

    counts = annotate_sbo_terms(model)

    assert counts["reactions"] == 0
    assert counts["metabolites"] == 0


# ---------------------------------------------------------------------------
# Integration test: enrich_model_annotations
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Unit tests: rename_to_bigg
# ---------------------------------------------------------------------------
def _real_model_for_rename():
    """Build a small real cobra.Model — needed because cobra's id setter
    enforces uniqueness through the model.dicts, which MagicMock can't fake.
    """
    import cobra

    m = cobra.Model("rename_test")
    # Two metabolites in c0 + one in e0, plus an exchange.
    for met_id in ("cpd00027_c0", "cpd00001_c0", "cpd00027_e0"):
        m.add_metabolites([cobra.Metabolite(met_id, compartment=met_id.split("_")[1])])
    # One internal reaction (rxn pattern) — gapseq-style with a bigg.reaction list
    r_internal = cobra.Reaction("rxn00001_c0")
    r_internal.add_metabolites(
        {
            m.metabolites.get_by_id("cpd00027_c0"): -1,
            m.metabolites.get_by_id("cpd00001_c0"): 1,
        }
    )
    r_internal.annotation["seed.reaction"] = "rxn00001"
    r_internal.annotation["bigg.reaction"] = ["PPA", "PPA_1pp"]
    m.add_reactions([r_internal])
    # One internal with NO existing bigg.reaction — TSV-only path
    r_internal_2 = cobra.Reaction("rxn00099_c0")
    r_internal_2.add_metabolites({m.metabolites.get_by_id("cpd00001_c0"): -1})
    r_internal_2.annotation["seed.reaction"] = "rxn00099"
    m.add_reactions([r_internal_2])
    # One exchange
    r_ex = cobra.Reaction("EX_cpd00027_e0")
    r_ex.add_metabolites({m.metabolites.get_by_id("cpd00027_e0"): -1})
    m.add_reactions([r_ex])
    return m


def test_rename_metabolite_with_compartment_translation():
    m = _real_model_for_rename()
    rename_to_bigg(m, seed_bigg={"cpd00027": "glc__D"}, seed_rxn_bigg={})
    assert "glc__D_c" in [met.id for met in m.metabolites]
    assert "cpd00027_c0" not in [met.id for met in m.metabolites]


def test_rename_exchange_pair():
    """The exchange reaction AND its boundary metabolite are renamed together."""
    m = _real_model_for_rename()
    rename_to_bigg(m, seed_bigg={"cpd00027": "glc__D"}, seed_rxn_bigg={})
    rxn_ids = {r.id for r in m.reactions}
    met_ids = {met.id for met in m.metabolites}
    assert "EX_glc__D_e" in rxn_ids
    assert "EX_cpd00027_e0" not in rxn_ids
    assert "glc__D_e" in met_ids


def test_rename_internal_reaction_via_tsv():
    """When no bigg.reaction annotation exists, fall back to the TSV."""
    m = _real_model_for_rename()
    rename_to_bigg(m, seed_bigg={}, seed_rxn_bigg={"rxn00099": "FOO_TSV_HIT"})
    rxn_ids = {r.id for r in m.reactions}
    assert "FOO_TSV_HIT" in rxn_ids
    assert "rxn00099_c0" not in rxn_ids


def test_rename_internal_reaction_prefers_existing_annotation():
    """Policy B: prefer first existing bigg.reaction over a TSV hit."""
    m = _real_model_for_rename()
    # The internal rxn00001 reaction has bigg.reaction = ['PPA', 'PPA_1pp']
    rename_to_bigg(m, seed_bigg={}, seed_rxn_bigg={"rxn00001": "WRONG_TSV_VALUE"})
    rxn_ids = {r.id for r in m.reactions}
    assert "PPA" in rxn_ids
    assert "WRONG_TSV_VALUE" not in rxn_ids


def test_rename_skips_unmapped_ids():
    m = _real_model_for_rename()
    counts = rename_to_bigg(m, seed_bigg={}, seed_rxn_bigg={})
    # No mappings provided -> only annotation-driven internal rename fires (rxn00001 -> PPA)
    assert counts["metabolites"] == 0
    assert counts["reactions_exchange"] == 0
    assert (
        counts["reactions_internal"] == 1
    )  # only the one with bigg.reaction annotation
    # Other entities keep SEED IDs
    assert "cpd00027_c0" in [m_.id for m_ in m.metabolites]
    assert "EX_cpd00027_e0" in [r.id for r in m.reactions]


def test_rename_skips_collision_with_existing_id():
    """Two cpd IDs that map to the same BiGG id collide; second skipped."""
    import cobra

    m = cobra.Model("c")
    m.add_metabolites(
        [
            cobra.Metabolite("cpd00027_c0", compartment="c0"),
            cobra.Metabolite("cpd99999_c0", compartment="c0"),
        ]
    )
    counts = rename_to_bigg(
        m,
        seed_bigg={"cpd00027": "glc__D", "cpd99999": "glc__D"},
        seed_rxn_bigg={},
    )
    met_ids = {met.id for met in m.metabolites}
    assert "glc__D_c" in met_ids
    assert ("cpd00027_c0" in met_ids) ^ (
        "cpd99999_c0" in met_ids
    )  # exactly one kept SEED
    assert counts["metabolites"] == 1
    assert counts["collisions_skipped"] == 1


def test_rename_is_idempotent():
    """A second rename pass is a no-op (regex doesn't match BiGG-form IDs)."""
    m = _real_model_for_rename()
    seed_bigg = {"cpd00027": "glc__D", "cpd00001": "h2o"}
    first = rename_to_bigg(m, seed_bigg=seed_bigg, seed_rxn_bigg={})
    second = rename_to_bigg(m, seed_bigg=seed_bigg, seed_rxn_bigg={})
    assert first["metabolites"] > 0
    assert second["metabolites"] == 0
    assert second["reactions_exchange"] == 0


def test_rename_preserves_seed_compound_annotation():
    """Pre-existing seed.compound / seed.reaction annotations survive renaming."""
    m = _real_model_for_rename()
    met = m.metabolites.get_by_id("cpd00027_c0")
    met.annotation["seed.compound"] = "cpd00027"  # gapseq native
    rename_to_bigg(m, seed_bigg={"cpd00027": "glc__D"}, seed_rxn_bigg={})
    renamed = m.metabolites.get_by_id("glc__D_c")
    assert renamed.annotation.get("seed.compound") == "cpd00027"
    # seed.reaction on the rxn00001 reaction stays after rename to PPA
    rxn = m.reactions.get_by_id("PPA")
    assert rxn.annotation.get("seed.reaction") == "rxn00001"


def test_rename_preserves_fba_growth_rate(tmp_path: Path):
    """Renaming must not change FBA biomass flux (cobrapy stoichiometry by reference)."""
    import cobra
    import cobra.io

    src = Path(__file__).parent / "data" / "tiny_model.xml"
    if not src.exists():
        pytest.skip("tests/data/tiny_model.xml fixture not present")
    m = cobra.io.read_sbml_model(str(src))
    try:
        biomass_before = m.slim_optimize()
    except Exception:
        pytest.skip(
            "tiny_model has no objective / not solvable; FBA invariance untested"
        )
    if biomass_before is None or biomass_before == 0:
        pytest.skip("tiny_model returns 0 / None biomass; FBA invariance untested")

    rename_to_bigg(m)
    biomass_after = m.slim_optimize()
    assert (
        abs(biomass_before - biomass_after) < 1e-6
    ), f"FBA invariance broken: {biomass_before} -> {biomass_after}"


def test_bare_ids_emit_miriam_uris_in_sbml_output(tmp_path: Path):
    """Pins the cobrapy convention: bare IDs in `annotation` are wrapped into
    https://identifiers.org/<ns>/<id> URIs by write_sbml_model. This is the
    reason annotate_* writes bare IDs (not URIs) into the dict — and the reason
    no explicit URI-emission step is needed in this module.
    """
    import cobra
    import cobra.io

    m = cobra.Model("uri_test")
    met = cobra.Metabolite("cpd00027_e0", compartment="e0")
    m.add_metabolites([met])
    rxn = cobra.Reaction("EX_cpd00027_e0")
    rxn.add_metabolites({met: -1})
    m.add_reactions([rxn])

    met.annotation["bigg.metabolite"] = "glc__D"
    met.annotation["seed.compound"] = "cpd00027"
    rxn.annotation["bigg.reaction"] = "EX_glc__D_e"
    rxn.annotation["seed.reaction"] = "rxn00001"

    out = tmp_path / "m.xml"
    cobra.io.write_sbml_model(m, str(out))
    text = out.read_text()

    for uri in (
        "https://identifiers.org/bigg.metabolite/glc__D",
        "https://identifiers.org/seed.compound/cpd00027",
        "https://identifiers.org/bigg.reaction/EX_glc__D_e",
        "https://identifiers.org/seed.reaction/rxn00001",
    ):
        assert uri in text, f"expected MIRIAM URI {uri!r} missing from SBML output"


def test_enrich_model_annotations(tmp_path: Path):
    """Full enrichment pipeline: genes + reactions + exchanges + internal metabolites."""
    in_file = tmp_path / "model.xml"
    out_file = tmp_path / "out" / "model.xml"

    genes = [_make_gene("AEP93129_1_1_308")]
    rxns = [_make_reaction("rxn1", {"ec-code": "3.6.1.1", "seed.reaction": "rxn00001"})]

    ex_met = MagicMock()
    ex_met.id = "cpd00027_e0"
    ex_met.annotation = {}
    ex_rxn = MagicMock()
    ex_rxn.id = "EX_cpd00027_e0"
    ex_rxn.annotation = {}
    ex_rxn.metabolites = {ex_met: -1.0}

    int_met = MagicMock()
    int_met.id = "cpd00082_c0"  # fructose, cytoplasmic
    int_met.annotation = {}

    mock_model = _make_model(genes=genes, reactions=rxns)
    mock_model.exchanges = [ex_rxn]
    mock_model.demands = []
    mock_model.sinks = []
    mock_model.metabolites = [ex_met, int_met]

    seed_rhea = {"rxn00001": "24576"}
    seed_bigg = {"cpd00027": "glc__D", "cpd00082": "fru"}

    def _mock_load(name):
        if "seed_to_rhea" in name:
            return seed_rhea
        if "seed_to_bigg" in name:
            return seed_bigg
        return {}

    with (
        patch(
            "parmmosahn.annotation.cobra.io.read_sbml_model", return_value=mock_model
        ),
        patch("parmmosahn.annotation.cobra.io.write_sbml_model") as mock_write,
        patch("parmmosahn.annotation._load_tsv_mapping", side_effect=_mock_load),
    ):
        summary = enrich_model_annotations(in_file, out_file)

    assert summary["genes_ncbiprotein"] == 1
    assert summary["reactions_brenda"] == 1
    assert summary["reactions_rhea"] == 1
    # exchange: bigg.metabolite on boundary met, bigg.reaction on reaction
    assert summary["exchanges_bigg"] == 1
    assert ex_met.annotation["bigg.metabolite"] == "glc__D"
    assert ex_rxn.annotation["bigg.reaction"] == "EX_glc__D_e"
    # intracellular: only int_met gets bigg.metabolite (ex_met already done)
    assert summary["metabolites_bigg"] == 1
    assert int_met.annotation["bigg.metabolite"] == "fru"
    # SBO terms: rxn1 (biochemical) + ex_met + int_met
    assert summary["sbo_reactions"] == 1
    assert summary["sbo_metabolites"] == 2
    mock_write.assert_called_once()


# ---------------------------------------------------------------------------
# CLI wiring test
# ---------------------------------------------------------------------------
def test_cli_annotate(tmp_path: Path):
    """'annotate' CLI command should wire through to enrich_models_directory."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "strain_A.xml").touch()
    out_dir = tmp_path / "annotated"
    runner = CliRunner()

    with patch("parmmosahn.annotation.enrich_models_directory") as mock_enrich:
        mock_enrich.return_value = Path(out_dir)
        result = runner.invoke(
            cli,
            ["annotate", "-M", str(models_dir), "-o", str(out_dir)],
        )

    assert result.exit_code == 0, result.output
    mock_enrich.assert_called_once()
