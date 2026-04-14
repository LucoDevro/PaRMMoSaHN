from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from parmmosahn.annotation import (
    _add_annotation,
    _get_annotation,
    annotate_genes,
    annotate_reactions,
    enrich_model_annotations,
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
# Integration test: enrich_model_annotations
# ---------------------------------------------------------------------------
def test_enrich_model_annotations(tmp_path: Path):
    """Full enrichment pipeline: genes + reactions."""
    in_file = tmp_path / "model.xml"
    out_file = tmp_path / "out" / "model.xml"

    genes = [_make_gene("AEP93129_1_1_308")]
    rxns = [_make_reaction("rxn1", {"ec-code": "3.6.1.1", "seed.reaction": "rxn00001"})]
    mock_model = _make_model(genes=genes, reactions=rxns)

    seed_rhea = {"rxn00001": "24576"}

    with (
        patch(
            "parmmosahn.annotation.cobra.io.read_sbml_model", return_value=mock_model
        ),
        patch("parmmosahn.annotation.cobra.io.write_sbml_model") as mock_write,
        patch(
            "parmmosahn.annotation._load_tsv_mapping",
            side_effect=lambda name: seed_rhea if "seed" in name else {},
        ),
    ):
        summary = enrich_model_annotations(in_file, out_file)

    assert summary["genes_ncbiprotein"] == 1
    assert summary["reactions_brenda"] == 1
    assert summary["reactions_rhea"] == 1
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
