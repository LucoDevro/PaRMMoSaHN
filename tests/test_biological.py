"""
Tier 3: Biological / Integration tests for PaRMMoSaHN.

These tests validate pre-computed pipeline artifacts against biological
invariants and detect known bugs. They require the ``test_run/`` directory
with pre-computed results from 4 Bacillus genomes.

Run with:  pytest -m biological
"""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl
import pytest

pytestmark = pytest.mark.biological

EXPECTED_STRAINS = [
    "GCA_000227485.1",
    "GCF_000009045.1",
    "GCF_000146565.1",
    "GCF_000209795.2",
]


def _require_strain(mapping: dict, strain: str):
    """Skip test if strain is not present in the pre-computed artifacts."""
    if strain not in mapping:
        pytest.skip(f"{strain}: not available in pre-computed artifacts")


# ---------------------------------------------------------------------------
# Biological range constants (derived from manual inspection of test_run)
# ---------------------------------------------------------------------------
MIN_REACTIONS = 2800
MAX_REACTIONS = 3200
MIN_METABOLITES = 2300
MAX_METABOLITES = 2700
MIN_GENES = 1100
MAX_GENES = 1400
EXPECTED_COMPARTMENTS = {"c0", "e0", "p0"}
MIN_EXCHANGE_REACTIONS = 180
EXPECTED_RXN_TABLE_COLUMNS = 19

# ===================================================================
# Section A: Pan-Results Validation (Steps 1-2 output)
# ===================================================================


class TestPanResults:
    """Validate pre-computed pangenome and pan-model artifacts."""

    def test_pangenome_fasta_valid(self, precomputed_pan_fasta: Path):
        text = precomputed_pan_fasta.read_text(encoding="utf-8")
        headers = [l for l in text.splitlines() if l.startswith(">")]
        ids = [h.split()[0][1:] for h in headers]

        assert len(ids) > 3000, f"Expected >3000 sequences, got {len(ids)}"
        assert len(ids) < 6000, f"Expected <6000 sequences, got {len(ids)}"
        assert len(ids) == len(set(ids)), "Duplicate sequence IDs in pan FASTA"

    def test_pangenome_h5_exists(self, test_run_dir: Path):
        h5 = test_run_dir / "00-pan_results" / "01-pangenome" / "pangenome.h5"
        assert h5.exists(), "pangenome.h5 not found"
        assert h5.stat().st_size > 1_000_000, "pangenome.h5 suspiciously small (<1MB)"

    def test_panmodel_reaction_table_valid(self, precomputed_pan_rxn_table: Path):
        df = pl.read_csv(
            precomputed_pan_rxn_table,
            separator="\t",
            comment_prefix="#",
            quote_char=None,
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        )
        assert (
            df.width == EXPECTED_RXN_TABLE_COLUMNS
        ), f"Expected {EXPECTED_RXN_TABLE_COLUMNS} columns, got {df.width}"
        assert len(df) > 50_000, f"Expected >50k rows, got {len(df)}"
        assert df["rxn"].null_count() == 0, "rxn column has nulls"
        assert "stitle" in df.columns, "stitle column missing"


# ===================================================================
# Section B: Strain Model Biological Validation (Step 3 output)
# ===================================================================


class TestStrainModelBiology:
    """Validate that strain SBML models are biologically plausible."""

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_strain_model_loads_as_valid_sbml(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        assert model is not None
        assert len(model.reactions) > 0

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_strain_model_size_in_range(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        assert (
            MIN_REACTIONS <= len(model.reactions) <= MAX_REACTIONS
        ), f"{strain}: {len(model.reactions)} reactions outside [{MIN_REACTIONS}, {MAX_REACTIONS}]"
        assert (
            MIN_METABOLITES <= len(model.metabolites) <= MAX_METABOLITES
        ), f"{strain}: {len(model.metabolites)} metabolites outside [{MIN_METABOLITES}, {MAX_METABOLITES}]"
        assert (
            MIN_GENES <= len(model.genes) <= MAX_GENES
        ), f"{strain}: {len(model.genes)} genes outside [{MIN_GENES}, {MAX_GENES}]"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_strain_model_has_three_compartments(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        assert set(model.compartments.keys()) == EXPECTED_COMPARTMENTS

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_strain_model_has_exchange_reactions(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        ex_rxns = [r for r in model.reactions if r.id.startswith("EX_")]
        assert (
            len(ex_rxns) >= MIN_EXCHANGE_REACTIONS
        ), f"{strain}: only {len(ex_rxns)} exchange reactions"
        r_ex_rxns = [r for r in model.reactions if r.id.startswith("R_EX_")]
        assert (
            len(r_ex_rxns) == 0
        ), f"{strain}: {len(r_ex_rxns)} reactions still have R_EX_ prefix"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_strain_model_has_biomass(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        bio_rxns = [r for r in model.reactions if r.id == "bio1"]
        assert len(bio_rxns) == 1, f"{strain}: expected 1 bio1 reaction"
        assert model.objective.expression is not None

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_strain_model_produces_biomass(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        solution = model.optimize()
        assert (
            solution.status == "optimal"
        ), f"{strain}: FBA status is {solution.status}, not optimal"
        assert (
            solution.objective_value > 0.1
        ), f"{strain}: biomass flux {solution.objective_value} too low"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_strain_model_essential_metabolites(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        met_ids = {m.id for m in model.metabolites}
        # cobrapy strips the M_ prefix from SBML metaid; actual IDs are cpd*_c0
        assert "cpd00001_c0" in met_ids, f"{strain}: water missing"
        assert "cpd00002_c0" in met_ids, f"{strain}: ATP missing"

    def test_strain_model_reaction_counts_vary(self, loaded_cobra_models):
        counts = [len(m.reactions) for m in loaded_cobra_models.values()]
        assert (
            len(set(counts)) > 1
        ), "All strains have identical reaction counts — models may be pan-model copies"


# ===================================================================
# Section C: Pipeline Artifact Integrity
# ===================================================================


class TestArtifactIntegrity:
    """Verify intermediate pipeline files exist and are well-formed."""

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_filtered_reactions_table_valid(
        self, precomputed_strain_dirs, precomputed_pan_rxn_table, strain
    ):
        _require_strain(precomputed_strain_dirs, strain)
        tbl = precomputed_strain_dirs[strain] / "filtered-Reactions.tbl"
        assert tbl.exists(), f"{strain}: filtered-Reactions.tbl missing"

        df = pl.read_csv(
            tbl,
            separator="\t",
            comment_prefix="#",
            quote_char=None,
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        )
        assert df.width == EXPECTED_RXN_TABLE_COLUMNS
        # Filtered table should have fewer rows than the pan-model
        df_pan = pl.read_csv(
            precomputed_pan_rxn_table,
            separator="\t",
            comment_prefix="#",
            quote_char=None,
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        )
        assert len(df) < len(
            df_pan
        ), f"{strain}: filtered table ({len(df)}) not smaller than pan table ({len(df_pan)})"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_diamond_matches_exist(self, strain_output_root, strain):
        tsv = strain_output_root / "matches" / f"{strain}.tsv"
        if not tsv.exists():
            pytest.skip(f"{strain}: DIAMOND match file not present")

        lines = tsv.read_text(encoding="utf-8").strip().splitlines()
        assert len(lines) > 1000, f"{strain}: only {len(lines)} DIAMOND hits"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_proteome_extracted(self, strain_output_root, strain):
        faa = strain_output_root / "proteomes" / f"{strain}.faa"
        if not faa.exists():
            pytest.skip(f"{strain}: proteome FAA not present")
        headers = [
            l for l in faa.read_text(encoding="utf-8").splitlines() if l.startswith(">")
        ]
        assert len(headers) > 1000, f"{strain}: only {len(headers)} proteins"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_filtered_reactions_stitle_has_real_ids(
        self, precomputed_strain_dirs, strain
    ):
        _require_strain(precomputed_strain_dirs, strain)
        tbl = precomputed_strain_dirs[strain] / "filtered-Reactions.tbl"
        df = pl.read_csv(
            tbl,
            separator="\t",
            comment_prefix="#",
            quote_char=None,
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        )
        stitle_values = df["stitle"].drop_nulls().to_list()
        assert len(stitle_values) > 0
        # At least 90% should look like protein accessions (NP_123.1, AEP89074.1, WP_123.1)
        accession_pat = re.compile(r"^[A-Za-z_]+\d+\.\d+$")
        matches = sum(1 for v in stitle_values if accession_pat.match(str(v)))
        ratio = matches / len(stitle_values)
        assert (
            ratio > 0.9
        ), f"{strain}: only {ratio:.0%} of stitle values look like protein accessions"


# ===================================================================
# Section D: Bug Detection Tests
# ===================================================================


class TestRegressions:
    """Regression tests for previously discovered bugs."""

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_bug_model_id_replaced(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        assert (
            model.id != "filtered_Reactions_tbl"
        ), f"{strain}: model ID still 'filtered_Reactions_tbl' (regex didn't match)"
        # Expected safe ID: dots/hyphens replaced with underscores, M_ prefix if starts with digit
        safe = re.sub(r"[^a-zA-Z0-9_]", "_", strain)
        if safe[0].isdigit():
            safe = f"M_{safe}"
        assert (
            model.id == safe
        ), f"{strain}: expected model ID '{safe}', got '{model.id}'"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_bug_gene_ids_are_biological(self, loaded_cobra_models, strain):
        _require_strain(loaded_cobra_models, strain)
        model = loaded_cobra_models[strain]
        gene_ids = [g.id for g in model.genes]
        # gapseq generates IDs like sequence_1015_1_308 (cobrapy strips the G_ prefix)
        sequence_pat = re.compile(r"^sequence_\d+(_\d+)*$")
        seq_count = sum(1 for gid in gene_ids if sequence_pat.match(gid))
        ratio = seq_count / len(gene_ids)
        assert (
            ratio < 0.1
        ), f"{strain}: {ratio:.0%} of gene IDs are sequence_N (expected <10%)"

    @pytest.mark.parametrize("strain", EXPECTED_STRAINS)
    def test_bug_gene_ids_match_reactions_table(
        self, loaded_cobra_models, precomputed_strain_dirs, strain
    ):
        _require_strain(loaded_cobra_models, strain)
        _require_strain(precomputed_strain_dirs, strain)
        model = loaded_cobra_models[strain]
        tbl = precomputed_strain_dirs[strain] / "filtered-Reactions.tbl"
        df = pl.read_csv(
            tbl,
            separator="\t",
            comment_prefix="#",
            quote_char=None,
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        )
        real_ids = set(df["stitle"].drop_nulls().unique().to_list())
        # Sanitize for SBML SId: dots → underscores (cobrapy strips the G_ prefix)
        sanitized_real = {re.sub(r"[^a-zA-Z0-9_]", "_", str(rid)) for rid in real_ids}
        model_gene_ids = {g.id for g in model.genes}
        # gapseq appends suffixes (e.g. NP_387883_1 becomes NP_387883_1_1_308)
        # so check if any model gene ID starts with a sanitized real ID
        matches = sum(
            1
            for gid in model_gene_ids
            if any(gid.startswith(rid) for rid in sanitized_real)
        )
        assert (
            matches > 0
        ), f"{strain}: no real protein IDs from filtered-Reactions.tbl found in model genes"


# ===================================================================
# Section E: Live Pipeline Execution Tests (slow, require tools)
# ===================================================================


@pytest.mark.slow
class TestLivePipeline:
    """Tests that execute actual pipeline steps on real data."""

    def test_live_filter_blast_table(
        self,
        precomputed_pan_rxn_table,
        strain_output_root,
        precomputed_strain_dirs,
        tmp_path,
    ):
        """Run filter_blast_table on real data and compare with pre-computed output."""
        from parmmosahn.models import filter_blast_table

        strain = "GCF_000009045.1"
        match_file = strain_output_root / "matches" / f"{strain}.tsv"
        if not match_file.exists():
            pytest.skip("DIAMOND match file not found")

        output = tmp_path / "filtered-Reactions.tbl"
        filter_blast_table(precomputed_pan_rxn_table, match_file, output)

        assert output.exists()
        df = pl.read_csv(
            output,
            separator="\t",
            quote_char=None,
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        )
        assert df.width == EXPECTED_RXN_TABLE_COLUMNS
        assert len(df) > 5000

        # Compare with pre-computed reference
        ref = precomputed_strain_dirs[strain] / "filtered-Reactions.tbl"
        df_ref = pl.read_csv(
            ref,
            separator="\t",
            quote_char=None,
            truncate_ragged_lines=True,
            infer_schema_length=10000,
        )
        assert len(df) == len(
            df_ref
        ), f"Row count mismatch: {len(df)} vs reference {len(df_ref)}"

    def test_live_derive_single_strain(
        self,
        require_pipeline_tools,
        precomputed_pan_fasta,
        precomputed_pan_rxn_table,
        strain_genomes_dir,
        medium_file,
        tmp_path,
    ):
        """Run derive_strain_models on a single genome end-to-end (~10 min)."""
        from parmmosahn.models import derive_strain_models

        # Copy just one genome to isolate the test
        single_genome_dir = tmp_path / "genomes"
        single_genome_dir.mkdir()
        src = strain_genomes_dir / "GCF_000009045.1.gbff"
        if not src.exists():
            pytest.skip("Reference genome GCF_000009045.1.gbff not found")
        import shutil

        shutil.copy(src, single_genome_dir / src.name)

        out = tmp_path / "output"
        models_dir = derive_strain_models(
            genomes_dir=single_genome_dir,
            pangenome_fasta=precomputed_pan_fasta,
            panmodel_rxn_table=precomputed_pan_rxn_table,
            output_dir=out,
            medium_file=medium_file,
            threads=4,
        )

        # Verify output structure
        strain = "GCF_000009045.1"
        assert (models_dir / strain / f"{strain}.xml").exists()
        assert (out / "matches" / f"{strain}.tsv").exists()
        assert (out / "proteomes" / f"{strain}.faa").exists()

        # Load and validate model
        import cobra as cb

        model = cb.io.read_sbml_model(str(models_dir / strain / f"{strain}.xml"))
        assert MIN_REACTIONS <= len(model.reactions) <= MAX_REACTIONS
        solution = model.optimize()
        assert solution.status == "optimal"
        assert solution.objective_value > 0.1
