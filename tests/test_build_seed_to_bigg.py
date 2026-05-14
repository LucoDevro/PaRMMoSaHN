"""Tests for scripts/build_seed_to_bigg.py — the MetaNetX -> TSV builder.

The build script is reachable as a module via sys.path manipulation; it is
not packaged into the wheel (it is a developer/maintainer utility), so the
test imports it from the repo's scripts/ directory directly.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "build_seed_to_bigg.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_seed_to_bigg", _SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def builder():
    return _load_module()


@pytest.fixture
def tiny_chem_xref(tmp_path):
    """Minimal chem_xref.tsv with header, two SEED rows, two BiGG rows."""
    p = tmp_path / "chem_xref.tsv"
    p.write_text(
        "#VERSION:   4.5\n"
        "#DATE:      2025/08/13\n"
        "#VERSION:   244 (other-resource version)\n"  # noise after first
        "#source\tID\tdescription\n"
        "seed.compound:cpd00001\tMNXM2\tH2O||water\n"
        "seed.compound:cpd00027\tMNXM41\tD-glucose\n"
        "seed.compound:cpd99999\tMNXM_UNMAPPED\tno BiGG link\n"
        "bigg.metabolite:h2o\tMNXM2\tH2O\n"
        "bigg.metabolite:glc__D\tMNXM41\tD-glucose\n"
        "bigg.metabolite:glc_alt\tMNXM41\tD-glucose alt entry\n"  # tests first-wins
    )
    return p


def test_parse_chem_xref_extracts_first_version_only(builder, tiny_chem_xref):
    _, _, version = builder.parse_chem_xref(tiny_chem_xref)
    assert version == "4.5"


def test_parse_chem_xref_collects_seed_and_bigg(builder, tiny_chem_xref):
    seed_to_mnx, mnx_to_bigg, _ = builder.parse_chem_xref(tiny_chem_xref)
    assert seed_to_mnx == {
        "cpd00001": "MNXM2",
        "cpd00027": "MNXM41",
        "cpd99999": "MNXM_UNMAPPED",
    }
    # First-wins for MNX -> BiGG (glc__D registered before glc_alt)
    assert mnx_to_bigg == {"MNXM2": "h2o", "MNXM41": "glc__D"}


def test_build_mapping_drops_seed_with_no_bigg_chain(builder):
    mapping = builder.build_mapping(
        seed_to_mnx={"cpd00001": "MNXM2", "cpd99999": "MNXM_X"},
        mnx_to_bigg={"MNXM2": "h2o"},
        overrides={},
    )
    assert mapping == {"cpd00001": "h2o"}


def test_build_mapping_overrides_win_over_chain(builder):
    """Overrides must replace chain-derived mappings, not be ignored."""
    mapping = builder.build_mapping(
        seed_to_mnx={"cpd00082": "MNXM1363862"},
        mnx_to_bigg={"MNXM1363862": "fruchose_wrong"},  # would chain to wrong BiGG
        overrides={"cpd00082": "fru"},
    )
    assert mapping["cpd00082"] == "fru"


def test_build_mapping_override_works_when_chain_misses(builder):
    """Overrides fire even when no chain exists for that SEED ID."""
    mapping = builder.build_mapping(
        seed_to_mnx={},
        mnx_to_bigg={},
        overrides={"cpd00611": "thr__D"},
    )
    assert mapping == {"cpd00611": "thr__D"}


def test_write_tsv_is_sorted_and_has_header(builder, tmp_path):
    out = tmp_path / "out.tsv"
    builder.write_tsv({"cpd00027": "glc__D", "cpd00001": "h2o"}, out)
    lines = out.read_text().splitlines()
    assert lines[0] == "seed_cpd_id\tbigg_met_id"
    assert lines[1] == "cpd00001\th2o"
    assert lines[2] == "cpd00027\tglc__D"


def test_write_tsv_deterministic(builder, tmp_path):
    """Same input produces the same bytes (no dict-iteration nondeterminism)."""
    out1 = tmp_path / "out1.tsv"
    out2 = tmp_path / "out2.tsv"
    mapping = {f"cpd{i:05d}": f"bigg_{i}" for i in range(50)}
    builder.write_tsv(mapping, out1)
    builder.write_tsv(mapping, out2)
    assert out1.read_bytes() == out2.read_bytes()


def test_overrides_dict_is_non_empty(builder):
    """Regression guard — the supervisor-verified overrides must not be lost."""
    assert "cpd00082" in builder.BIGG_OVERRIDES  # D-Fructose
    assert builder.BIGG_OVERRIDES["cpd00082"] == "fru"
    assert "cpd00611" in builder.BIGG_OVERRIDES  # D-Threonine
    assert len(builder.BIGG_OVERRIDES) >= 14
