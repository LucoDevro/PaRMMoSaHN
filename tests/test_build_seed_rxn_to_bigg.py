"""Tests for scripts/build_seed_rxn_to_bigg.py — MetaNetX reaction xref -> TSV."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_SCRIPT = (
    Path(__file__).resolve().parent.parent / "scripts" / "build_seed_rxn_to_bigg.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location("build_seed_rxn_to_bigg", _SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def builder():
    return _load_module()


@pytest.fixture
def tiny_reac_xref(tmp_path):
    p = tmp_path / "reac_xref.tsv"
    p.write_text(
        "#VERSION:   4.5\n"
        "#DATE:      2025/08/13\n"
        "#VERSION:   1.6.0 (other-resource)\n"
        "#source\tID\tdescription\n"
        "seed.reaction:rxn00001\tMNXR1\tpyrophosphatase\n"
        "seed.reaction:rxn00007\tMNXR7\ttrehalose hydrolysis\n"
        "seed.reaction:rxn99999\tMNXR_X\tunmapped\n"
        "bigg.reaction:IDP_h\tMNXR1\tinorganic diphosphatase (chloroplast)\n"
        "bigg.reaction:PPA\tMNXR1\tpyrophosphatase\n"  # tests first-wins
        "bigg.reaction:R_TREH\tMNXR7\ttrehalase\n"  # tests R_ stripping
    )
    return p


def test_parse_extracts_version_and_dicts(builder, tiny_reac_xref):
    seed, mnxr, version = builder.parse_reac_xref(tiny_reac_xref)
    assert version == "4.5"
    assert seed == {
        "rxn00001": "MNXR1",
        "rxn00007": "MNXR7",
        "rxn99999": "MNXR_X",
    }


def test_r_prefix_is_stripped_on_bigg_ids(builder, tiny_reac_xref):
    """SBML SId-safe `R_` prefix is removed from harvested BiGG reaction IDs."""
    _, mnxr_to_bigg, _ = builder.parse_reac_xref(tiny_reac_xref)
    assert mnxr_to_bigg["MNXR7"] == "TREH"  # not "R_TREH"


def test_first_wins_for_multiple_bigg_per_mnxr(builder, tiny_reac_xref):
    """When two BiGG IDs map to the same MNXR, the first one is kept."""
    _, mnxr_to_bigg, _ = builder.parse_reac_xref(tiny_reac_xref)
    assert mnxr_to_bigg["MNXR1"] == "IDP_h"  # first in source order


def test_build_mapping_chains_seed_to_bigg(builder, tiny_reac_xref):
    seed, mnxr_to_bigg, _ = builder.parse_reac_xref(tiny_reac_xref)
    mapping = builder.build_mapping(seed, mnxr_to_bigg, overrides={})
    assert mapping == {
        "rxn00001": "IDP_h",
        "rxn00007": "TREH",
    }
    assert "rxn99999" not in mapping


def test_overrides_take_precedence(builder):
    mapping = builder.build_mapping(
        seed_to_mnxr={"rxn00001": "MNXR1"},
        mnxr_to_bigg={"MNXR1": "wrong_default"},
        overrides={"rxn00001": "PPA"},
    )
    assert mapping["rxn00001"] == "PPA"


def test_write_tsv_sorted_with_header(builder, tmp_path):
    out = tmp_path / "out.tsv"
    builder.write_tsv({"rxn00007": "TREH", "rxn00001": "IDP_h"}, out)
    lines = out.read_text().splitlines()
    assert lines[0] == "seed_rxn_id\tbigg_rxn_id"
    assert lines[1] == "rxn00001\tIDP_h"
    assert lines[2] == "rxn00007\tTREH"


def test_overrides_dict_exists(builder):
    """The overrides dict is present (may be empty initially)."""
    assert hasattr(builder, "BIGG_RXN_OVERRIDES")
    assert isinstance(builder.BIGG_RXN_OVERRIDES, dict)
