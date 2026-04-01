"""Tests for the optional analysis extensions (analysis.py) and their CLI wiring."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from parmmosahn.cli import cli

# ---------------------------------------------------------------------------
# Helpers: build minimal cobra-like model mocks
# ---------------------------------------------------------------------------


def _make_mock_model(name: str, reaction_ids: list[str], growth: float = 0.5):
    """Return a MagicMock that quacks like a cobra.Model."""
    model = MagicMock()
    model.id = name

    # Reactions — build mock objects, expose as an iterable container
    rxns = []
    for rid in reaction_ids:
        r = MagicMock()
        r.id = rid
        rxns.append(r)

    rxn_container = MagicMock()
    rxn_container.__iter__ = MagicMock(return_value=iter(rxns))
    model.reactions = rxn_container

    # Exchange reactions (prefix EX_)
    exchanges = [r for r in rxns if r.id.startswith("EX_")]
    model.exchanges = exchanges

    # Medium — one entry per exchange
    model.medium = {r.id: 10.0 for r in exchanges}

    # optimize() returns a solution
    sol = MagicMock()
    sol.status = "optimal"
    sol.objective_value = growth
    sol.fluxes = MagicMock()
    sol.fluxes.get = lambda rid, default=0.0: (-5.0 if rid.startswith("EX_") else 0.0)
    model.optimize = MagicMock(return_value=sol)

    return model


# ---------------------------------------------------------------------------
# Unit tests — fba_summary
# ---------------------------------------------------------------------------


@patch("parmmosahn.analysis._load_models")
def test_fba_summary_writes_tsv(mock_load, tmp_path: Path):
    """fba_summary should produce a TSV with one header + one row per model."""
    from parmmosahn.analysis import fba_summary

    m1 = _make_mock_model("strain_A", ["RXN_1", "EX_glc", "EX_o2"], growth=0.42)
    m2 = _make_mock_model("strain_B", ["RXN_1", "RXN_2", "EX_glc"], growth=0.38)
    mock_load.return_value = [("strain_A", m1), ("strain_B", m2)]

    out = tmp_path / "fba.tsv"
    result = fba_summary(tmp_path, out)

    assert result == out
    lines = out.read_text().strip().splitlines()
    assert lines[0].startswith("strain\t")
    assert len(lines) == 3  # header + 2 strains
    assert "strain_A" in lines[1]
    assert "0.42" in lines[1]


@patch("parmmosahn.analysis._load_models")
def test_fba_summary_infeasible(mock_load, tmp_path: Path):
    """Infeasible model should show status but no active exchanges."""
    from parmmosahn.analysis import fba_summary

    m1 = _make_mock_model("bad", ["EX_glc"], growth=0.0)
    sol = m1.optimize.return_value
    sol.status = "infeasible"
    sol.objective_value = None
    mock_load.return_value = [("bad", m1)]

    out = tmp_path / "fba.tsv"
    fba_summary(tmp_path, out)

    content = out.read_text()
    assert "infeasible" in content
    assert "NA" in content


# ---------------------------------------------------------------------------
# Unit tests — auxotrophy_screen
# ---------------------------------------------------------------------------


@patch("parmmosahn.analysis._load_models")
def test_auxotrophy_screen_matrix(mock_load, tmp_path: Path):
    """auxotrophy_screen should write a component x strain matrix + baseline row."""
    from parmmosahn.analysis import auxotrophy_screen

    m1 = _make_mock_model("s1", ["EX_glc", "EX_nh4"], growth=0.5)
    m2 = _make_mock_model("s2", ["EX_glc", "EX_nh4"], growth=0.6)
    mock_load.return_value = [("s1", m1), ("s2", m2)]

    out = tmp_path / "auxo.tsv"
    auxotrophy_screen(tmp_path, out)

    lines = out.read_text().strip().splitlines()
    header = lines[0].split("\t")
    assert header[0] == "component"
    assert "s1" in header and "s2" in header

    # Should have: header + N components + _baseline row
    # Both models share EX_glc + EX_nh4 → 2 components + 1 baseline + 1 header = 4
    assert len(lines) == 4
    assert lines[-1].startswith("_baseline")


@patch("parmmosahn.analysis._load_models")
def test_auxotrophy_detects_essential(mock_load, tmp_path: Path):
    """When removing a component kills growth, it should be logged as auxotrophic."""
    from parmmosahn.analysis import auxotrophy_screen

    model = _make_mock_model("s1", ["EX_glc", "EX_nh4"], growth=0.5)

    # Make knockouts return different growth values
    call_count = {"n": 0}
    original_optimize = model.optimize

    def side_effect_optimize():
        sol = MagicMock()
        sol.status = "optimal"
        # First call is baseline (0.5), then EX_glc knockout → 0, EX_nh4 knockout → 0.4
        values = [0.5, 0.0, 0.4]
        sol.objective_value = values[min(call_count["n"], len(values) - 1)]
        sol.fluxes = MagicMock()
        sol.fluxes.get = lambda rid, d=0.0: 0.0
        call_count["n"] += 1
        return sol

    model.optimize = side_effect_optimize
    mock_load.return_value = [("s1", model)]

    out = tmp_path / "auxo.tsv"
    auxotrophy_screen(tmp_path, out)

    content = out.read_text()
    assert "0.000000" in content  # EX_glc knockout → zero growth


# ---------------------------------------------------------------------------
# Unit tests — reaction_heatmap
# ---------------------------------------------------------------------------


@patch("parmmosahn.analysis._load_models")
def test_reaction_heatmap_matrix(mock_load, tmp_path: Path):
    """reaction_heatmap should write a binary matrix of reactions x strains."""
    from parmmosahn.analysis import reaction_heatmap

    m1 = _make_mock_model("s1", ["RXN_A", "RXN_B", "RXN_C"])
    m2 = _make_mock_model("s2", ["RXN_A", "RXN_C", "RXN_D"])
    mock_load.return_value = [("s1", m1), ("s2", m2)]

    out = tmp_path / "heatmap.tsv"
    reaction_heatmap(tmp_path, out)

    lines = out.read_text().strip().splitlines()
    header = lines[0].split("\t")
    assert header == ["reaction", "s1", "s2"]

    # 4 unique reactions
    assert len(lines) == 5  # header + 4 reactions

    # Check RXN_B: present in s1, absent in s2
    rxn_b_line = [l for l in lines if l.startswith("RXN_B")][0]
    parts = rxn_b_line.split("\t")
    assert parts[1] == "1" and parts[2] == "0"

    # Check RXN_D: absent in s1, present in s2
    rxn_d_line = [l for l in lines if l.startswith("RXN_D")][0]
    parts = rxn_d_line.split("\t")
    assert parts[1] == "0" and parts[2] == "1"


@patch("parmmosahn.analysis._load_models")
def test_reaction_heatmap_no_plot_without_flag(mock_load, tmp_path: Path):
    """plot=False should not produce a PNG."""
    from parmmosahn.analysis import reaction_heatmap

    m1 = _make_mock_model("s1", ["RXN_A"])
    mock_load.return_value = [("s1", m1)]

    out = tmp_path / "heatmap.tsv"
    reaction_heatmap(tmp_path, out, plot=False)

    assert out.exists()
    assert not out.with_suffix(".png").exists()


# ---------------------------------------------------------------------------
# CLI wiring tests
# ---------------------------------------------------------------------------


def test_cli_analyze_fba(tmp_path: Path):
    """'analyze fba' should wire through to fba_summary."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "s1.xml").touch()
    out = tmp_path / "fba.tsv"
    runner = CliRunner()

    with patch("parmmosahn.analysis.fba_summary") as mock_fba:
        mock_fba.return_value = out
        result = runner.invoke(
            cli,
            ["analyze", "fba", "-M", str(models_dir), "-o", str(out)],
        )

    assert result.exit_code == 0, result.output
    mock_fba.assert_called_once()


def test_cli_analyze_auxotrophy(tmp_path: Path):
    """'analyze auxotrophy' should wire through to auxotrophy_screen."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "s1.xml").touch()
    out = tmp_path / "auxo.tsv"
    runner = CliRunner()

    with patch("parmmosahn.analysis.auxotrophy_screen") as mock_auxo:
        mock_auxo.return_value = out
        result = runner.invoke(
            cli,
            ["analyze", "auxotrophy", "-M", str(models_dir), "-o", str(out)],
        )

    assert result.exit_code == 0, result.output
    mock_auxo.assert_called_once()


def test_cli_analyze_heatmap(tmp_path: Path):
    """'analyze heatmap' should wire through to reaction_heatmap."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "s1.xml").touch()
    out = tmp_path / "heat.tsv"
    runner = CliRunner()

    with patch("parmmosahn.analysis.reaction_heatmap") as mock_heat:
        mock_heat.return_value = out
        result = runner.invoke(
            cli,
            ["analyze", "heatmap", "-M", str(models_dir), "-o", str(out)],
        )

    assert result.exit_code == 0, result.output
    mock_heat.assert_called_once()


def test_cli_analyze_heatmap_with_plot(tmp_path: Path):
    """'analyze heatmap --plot' should pass plot=True."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "s1.xml").touch()
    out = tmp_path / "heat.tsv"
    runner = CliRunner()

    with patch("parmmosahn.analysis.reaction_heatmap") as mock_heat:
        mock_heat.return_value = out
        result = runner.invoke(
            cli,
            ["analyze", "heatmap", "-M", str(models_dir), "-o", str(out), "--plot"],
        )

    assert result.exit_code == 0, result.output
    _, kwargs = mock_heat.call_args
    assert kwargs["plot"] is True


def test_cli_analyze_auxotrophy_with_threshold(tmp_path: Path):
    """'analyze auxotrophy --threshold' should forward custom threshold."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "s1.xml").touch()
    out = tmp_path / "auxo.tsv"
    runner = CliRunner()

    with patch("parmmosahn.analysis.auxotrophy_screen") as mock_auxo:
        mock_auxo.return_value = out
        result = runner.invoke(
            cli,
            [
                "analyze",
                "auxotrophy",
                "-M",
                str(models_dir),
                "-o",
                str(out),
                "--threshold",
                "0.01",
            ],
        )

    assert result.exit_code == 0, result.output
    _, kwargs = mock_auxo.call_args
    assert kwargs["threshold"] == 0.01
