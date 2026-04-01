from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from parmmosahn.cli import cli

# ---------------------------------------------------------------------------
# Patch targets
#
# Functions imported at module level in cli.py → patch as parmmosahn.cli.*
#   - derive_strain_models, build_pan_model, build_pangenome
#
# Functions imported lazily inside command bodies → patch at source module
#   - run_pan, run_downstream  → parmmosahn.pipeline.*
#   - generate_curation_template → parmmosahn.curation.*
#   - check_environment → parmmosahn.utils.*
#   - load_config → parmmosahn.utils.*
# ---------------------------------------------------------------------------


def test_cli_pan(tmp_path: Path):
    """'pan' command should wire arguments through to run_pan."""
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    runner = CliRunner()

    with patch("parmmosahn.pipeline.run_pan") as mock_run_pan:
        result = runner.invoke(
            cli,
            [
                "pan",
                "-g",
                str(genomes_dir),
                "-l",
                "clostridiales",
                "-t",
                "4",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_run_pan.assert_called_once_with(
        genomes_dir,
        Path("00-pan_results"),
        "clostridiales",
        threads=4,
        add_unique=False,
        soft_core_override=None,
    )


def test_cli_pan_missing_required_args():
    """'pan' should fail when --genomes or --label are missing."""
    runner = CliRunner()

    result = runner.invoke(cli, ["pan", "-l", "test"])
    assert result.exit_code != 0
    assert "Missing option" in result.output

    result = runner.invoke(cli, ["pan", "-g", "/tmp"])
    assert result.exit_code != 0
    assert "Missing option" in result.output


def test_cli_strains(tmp_path: Path):
    """'strains' should forward all options including new performance knobs."""
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    fasta = tmp_path / "pan.faa"
    fasta.touch()
    rxns = tmp_path / "pan.tbl"
    rxns.touch()
    medium = tmp_path / "medium.csv"
    medium.touch()
    runner = CliRunner()

    with patch("parmmosahn.cli.derive_strain_models") as mock_derive:
        result = runner.invoke(
            cli,
            [
                "strains",
                "-g",
                str(genomes_dir),
                "-f",
                str(fasta),
                "-r",
                str(rxns),
                "-m",
                str(medium),
                "-t",
                "8",
                "--diamond-bits",
                "200",
                "--gapseq-bits",
                "180",
                "--biomass",
                "neg",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_derive.assert_called_once_with(
        genomes_dir=genomes_dir,
        pangenome_fasta=fasta,
        panmodel_rxn_table=rxns,
        output_dir=Path("03-strain_models"),
        medium_file=medium,
        threads=8,
        diamond_bits=200,
        gapseq_bits=180,
        biomass="neg",
    )


def test_cli_strains_default_performance_knobs(tmp_path: Path):
    """Performance knob defaults should be 150/150/pos."""
    genomes_dir = tmp_path / "g"
    genomes_dir.mkdir()
    fasta = tmp_path / "f.faa"
    fasta.touch()
    rxns = tmp_path / "r.tbl"
    rxns.touch()
    medium = tmp_path / "m.csv"
    medium.touch()
    runner = CliRunner()

    with patch("parmmosahn.cli.derive_strain_models") as mock_derive:
        runner.invoke(
            cli,
            [
                "strains",
                "-g",
                str(genomes_dir),
                "-f",
                str(fasta),
                "-r",
                str(rxns),
                "-m",
                str(medium),
            ],
        )

    _, kwargs = mock_derive.call_args
    assert kwargs["diamond_bits"] == 150
    assert kwargs["gapseq_bits"] == 150
    assert kwargs["biomass"] == "pos"


def test_cli_project_normal_mode(tmp_path: Path):
    """'project' in normal mode calls derive_strain_models then run_downstream."""
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    pangenome = tmp_path / "pan.faa"
    pangenome.touch()
    pan_rxns = tmp_path / "pan.tbl"
    pan_rxns.touch()
    medium = tmp_path / "medium.csv"
    medium.touch()
    out = tmp_path / "out"
    runner = CliRunner()

    with (
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.pipeline.run_downstream") as mock_downstream,
    ):
        mock_derive.return_value = out / "03-strains" / "models"
        result = runner.invoke(
            cli,
            [
                "project",
                "-g",
                str(genomes),
                "-f",
                str(pangenome),
                "-r",
                str(pan_rxns),
                "-m",
                str(medium),
                "-o",
                str(out),
            ],
        )

    assert result.exit_code == 0, result.output
    mock_derive.assert_called_once()
    mock_downstream.assert_called_once()


def test_cli_project_performance_knobs_forwarded(tmp_path: Path):
    """'project' should forward --diamond-bits/--gapseq-bits/--biomass."""
    genomes = tmp_path / "g"
    genomes.mkdir()
    fasta = tmp_path / "f.faa"
    fasta.touch()
    rxns = tmp_path / "r.tbl"
    rxns.touch()
    medium = tmp_path / "m.csv"
    medium.touch()
    out = tmp_path / "out"
    runner = CliRunner()

    with (
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.pipeline.run_downstream"),
    ):
        mock_derive.return_value = out / "03-strains" / "models"
        runner.invoke(
            cli,
            [
                "project",
                "-g",
                str(genomes),
                "-f",
                str(fasta),
                "-r",
                str(rxns),
                "-m",
                str(medium),
                "-o",
                str(out),
                "--diamond-bits",
                "120",
                "--gapseq-bits",
                "90",
                "--biomass",
                "neg",
            ],
        )

    _, kwargs = mock_derive.call_args
    assert kwargs["diamond_bits"] == 120
    assert kwargs["gapseq_bits"] == 90
    assert kwargs["biomass"] == "neg"


def test_cli_project_resume_mode(tmp_path: Path):
    """'project --resume' should call run_downstream with the existing curation DB."""
    out = tmp_path / "out"
    drafts = out / "03-strains" / "models"
    drafts.mkdir(parents=True)
    curation_db = out / "curation_template.xlsx"
    curation_db.touch()
    runner = CliRunner()

    with patch("parmmosahn.pipeline.run_downstream") as mock_downstream:
        result = runner.invoke(
            cli,
            ["project", "-o", str(out), "--resume"],
        )

    assert result.exit_code == 0, result.output
    mock_downstream.assert_called_once()
    call_kwargs = mock_downstream.call_args.kwargs
    assert call_kwargs["curation_db"] == curation_db


def test_cli_project_resume_missing_drafts(tmp_path: Path):
    """'project --resume' should fail gracefully when draft dir doesn't exist."""
    out = tmp_path / "out"
    out.mkdir()
    runner = CliRunner()

    result = runner.invoke(cli, ["project", "-o", str(out), "--resume"])
    assert result.exit_code != 0


def test_cli_project_pause_for_curation(tmp_path: Path):
    """--pause_for_curation should stop after generating template and echo abs path."""
    genomes = tmp_path / "g"
    genomes.mkdir()
    fasta = tmp_path / "f.faa"
    fasta.touch()
    rxns = tmp_path / "r.tbl"
    rxns.touch()
    medium = tmp_path / "m.csv"
    medium.touch()
    out = tmp_path / "out"
    runner = CliRunner()

    with (
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.curation.generate_curation_template"),
        patch("parmmosahn.pipeline.run_downstream") as mock_downstream,
    ):
        mock_derive.return_value = out / "03-strains" / "models"
        result = runner.invoke(
            cli,
            [
                "project",
                "-g",
                str(genomes),
                "-f",
                str(fasta),
                "-r",
                str(rxns),
                "-m",
                str(medium),
                "-o",
                str(out),
                "--pause_for_curation",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_downstream.assert_not_called()
    assert "curation_template.xlsx" in result.output
    assert str(out.resolve()) in result.output


def test_cli_evaluate(tmp_path: Path):
    """'evaluate' should run memote on each XML model."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "strain_A.xml").touch()
    (models_dir / "strain_B.xml").touch()
    out_dir = tmp_path / "reports"
    runner = CliRunner()

    with patch("parmmosahn.curation.run_memote_snapshot") as mock_memote:
        result = runner.invoke(
            cli,
            ["evaluate", "-M", str(models_dir), "-o", str(out_dir)],
        )

    assert result.exit_code == 0, result.output
    assert mock_memote.call_count == 2
    assert "Evaluating 2 model(s)" in result.output


def test_cli_evaluate_with_curation_template(tmp_path: Path):
    """'evaluate --curation-template' should also generate the Excel file."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "strain_A.xml").touch()
    out_dir = tmp_path / "reports"
    runner = CliRunner()

    with (
        patch("parmmosahn.curation.run_memote_snapshot"),
        patch("parmmosahn.curation.extract_memote_problems") as mock_extract,
    ):
        result = runner.invoke(
            cli,
            [
                "evaluate",
                "-M",
                str(models_dir),
                "-o",
                str(out_dir),
                "--curation-template",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_extract.assert_called_once()
    assert "Curation template saved" in result.output


def test_cli_doctor():
    """'doctor' should report OK/MISSING and print a summary."""
    runner = CliRunner()

    fake_results = {
        "ppanggolin": {
            "found": True,
            "version": "ppanggolin 2.2.6",
            "path": "/usr/bin/ppanggolin",
        },
        "gapseq": {
            "found": True,
            "version": "gapseq version 1.4",
            "path": "/usr/bin/gapseq",
        },
        "diamond": {"found": False, "version": None, "path": None},
        "docker": {"found": True, "version": "Docker 24.0", "path": "/usr/bin/docker"},
        "podman": {"found": False, "version": None, "path": None},
        "apptainer": {"found": False, "version": None, "path": None},
    }

    with patch("parmmosahn.utils.check_environment", return_value=fake_results):
        result = runner.invoke(cli, ["doctor"])

    assert result.exit_code == 0
    assert "[OK]" in result.output
    assert "[MISSING]" in result.output
    assert "diamond" in result.output
    assert "WARNING" in result.output  # diamond missing → warning


def test_cli_init_config(tmp_path: Path):
    """'init-config' should write a YAML template to the specified path."""
    runner = CliRunner()
    out_file = tmp_path / "my_config.yml"

    result = runner.invoke(cli, ["init-config", "-o", str(out_file)])

    assert result.exit_code == 0, result.output
    assert out_file.exists()
    content = out_file.read_text()
    assert "threads" in content
    assert "diamond_bits" in content
    assert "biomass" in content


def test_cli_project_config_file(tmp_path: Path):
    """'project --config' should load defaults from the YAML file."""
    pytest_importorskip_yaml()

    genomes = tmp_path / "g"
    genomes.mkdir()
    fasta = tmp_path / "f.faa"
    fasta.touch()
    rxns = tmp_path / "r.tbl"
    rxns.touch()
    medium = tmp_path / "m.csv"
    medium.touch()
    out = tmp_path / "out"

    config_file = tmp_path / "config.yml"
    config_file.write_text(
        f"genomes: {genomes}\nfasta: {fasta}\npanmodel_rxns: {rxns}\n"
        f"medium: {medium}\nthreads: 6\ndiamond_bits: 111\n",
        encoding="utf-8",
    )

    runner = CliRunner()

    with (
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.pipeline.run_downstream"),
    ):
        mock_derive.return_value = out / "03-strains" / "models"
        result = runner.invoke(
            cli,
            ["project", "-o", str(out), "--config", str(config_file)],
        )

    assert result.exit_code == 0, result.output
    mock_derive.assert_called_once()
    _, kwargs = mock_derive.call_args
    assert kwargs["threads"] == 6
    assert kwargs["diamond_bits"] == 111


def pytest_importorskip_yaml():
    """Skip test if pyyaml is not installed."""
    import pytest

    pytest.importorskip(
        "yaml", reason="pyyaml not installed (pip install parmmosahn[config])"
    )
