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
        Path("."),
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
        parallel_strains=None,
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
        mock_derive.return_value = out / "03-strain_models" / "models"
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
        mock_derive.return_value = out / "03-strain_models" / "models"
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
    drafts = out / "03-strain_models" / "models"
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
        patch("parmmosahn.curation.generate_curation_template") as mock_gen,
        patch("parmmosahn.pipeline.run_downstream") as mock_downstream,
    ):
        mock_derive.return_value = out / "03-strain_models" / "models"
        mock_gen.return_value = out / "curation_template.xlsx"
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
    assert "curation_template" in result.output
    assert str(out.resolve()) in result.output


def test_cli_project_pause_for_curation_tty_continues_inline(tmp_path: Path):
    """In a TTY, --pause_for_curation prompts then continues to downstream."""
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

    template_path = out / "curation_template.xlsx"

    with (
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.curation.generate_curation_template") as mock_gen,
        patch("parmmosahn.pipeline.run_downstream") as mock_downstream,
        patch("parmmosahn.cli._stdin_is_tty", return_value=True),
        # click.pause reads a key from stdin — short-circuit it to a no-op.
        patch("click.pause") as mock_pause,
    ):
        mock_derive.return_value = out / "03-strain_models" / "models"
        mock_gen.return_value = template_path
        mock_downstream.return_value = {}
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
    mock_pause.assert_called_once()
    # After Enter, run_downstream should fire with the generated template as
    # the curation_db argument, not None.
    mock_downstream.assert_called_once()
    _, kwargs = mock_downstream.call_args
    assert kwargs["curation_db"] == template_path


def test_cli_project_pause_for_curation_tty_abort(tmp_path: Path):
    """Ctrl-C during the TTY pause prints a resume hint and exits cleanly."""
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

    import click as _click

    with (
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.curation.generate_curation_template") as mock_gen,
        patch("parmmosahn.pipeline.run_downstream") as mock_downstream,
        patch("parmmosahn.cli._stdin_is_tty", return_value=True),
        patch("click.pause", side_effect=_click.Abort()),
    ):
        mock_derive.return_value = out / "03-strain_models" / "models"
        mock_gen.return_value = out / "curation_template.xlsx"
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
    assert "Resume with:" in result.output or "Aborted" in result.output


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
    # A3: host resources block must appear so users see RAM before launching
    # memory-hungry steps (memote).
    assert "Host resources:" in result.output
    # #7: bundled data files block lists the runtime TSV lookups.
    assert "Bundled data files" in result.output
    assert "seed_to_bigg.tsv" in result.output


def test_echo_and_log_writes_to_file_handler(tmp_path: Path):
    """_echo_and_log must mirror its message into any FileHandler on the root
    logger (the mechanism run/project use to capture status into run.log)."""
    import logging

    from parmmosahn.cli import _echo_and_log

    log_path = tmp_path / "run.log"
    handler = logging.FileHandler(str(log_path), mode="w")
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(message)s"))
    root = logging.getLogger()
    root.addHandler(handler)
    try:
        _echo_and_log("hello from a pipeline step")
        handler.flush()
        assert "hello from a pipeline step" in log_path.read_text()
    finally:
        root.removeHandler(handler)
        handler.close()


def test_echo_and_log_no_file_handler_is_a_noop_beyond_echo():
    """Without a FileHandler the helper just behaves like click.echo —
    no exception, no implicit logger.info on the root StreamHandler."""
    from parmmosahn.cli import _echo_and_log

    # No FileHandler installed; this should not raise.
    _echo_and_log("standalone command output")


def test_cli_doctor_json_output():
    """doctor --json emits a single parseable JSON object on stdout."""
    import json as _json

    runner = CliRunner()
    fake_env = {
        "ppanggolin": {
            "found": True,
            "version": "2.2.6",
            "path": "/usr/bin/ppanggolin",
        },
        "gapseq": {"found": True, "version": "1.4.2", "path": "/usr/bin/gapseq"},
        "diamond": {"found": True, "version": "2.1.8", "path": "/usr/bin/diamond"},
        "docker": {"found": True, "version": "24.0", "path": "/usr/bin/docker"},
        "podman": {"found": False, "version": None, "path": None},
        "apptainer": {"found": False, "version": None, "path": None},
    }
    fake_data = {
        "mnxr_to_rhea.tsv": {
            "found": True,
            "size_bytes": 100,
            "rows": 10,
            "path": "/x",
        },
        "seed_to_rhea.tsv": {
            "found": True,
            "size_bytes": 100,
            "rows": 10,
            "path": "/x",
        },
        "seed_to_bigg.tsv": {
            "found": True,
            "size_bytes": 200,
            "rows": 20,
            "path": "/x",
        },
    }

    fake_resources = {
        "platform": "Linux-x86_64",
        "python_version": "3.12.0",
        "cpu_count": 16,
        "total_ram_gb": 32.0,
        "available_ram_gb": 24.0,
    }

    with (
        patch("parmmosahn.utils.check_environment", return_value=fake_env),
        patch("parmmosahn.utils.check_bundled_data", return_value=fake_data),
        patch("parmmosahn.utils.get_system_resources", return_value=fake_resources),
    ):
        result = runner.invoke(cli, ["doctor", "--json"])

    assert result.exit_code == 0, result.output
    # In JSON mode, NONE of the human-readable headers should appear.
    assert "PaRMMoSaHN Environment Check" not in result.output
    assert "Required tools:" not in result.output

    # The output must be a single valid JSON object covering all blocks.
    payload = _json.loads(result.output)
    assert set(payload["required_tools"]) == {"ppanggolin", "gapseq", "diamond"}
    assert set(payload["container_engines"]) == {"docker", "podman", "apptainer"}
    assert set(payload["bundled_data"]) == {
        "mnxr_to_rhea.tsv",
        "seed_to_rhea.tsv",
        "seed_to_bigg.tsv",
    }
    assert payload["checks"]["all_required_tools_ok"] is True
    assert payload["checks"]["any_container_engine_ok"] is True
    assert payload["checks"]["all_bundled_data_ok"] is True
    assert payload["warnings"] == []  # nothing wrong → no warnings
    assert "parmmosahn_version" in payload
    assert "host_resources" in payload


def test_cli_doctor_json_surfaces_warnings_as_list():
    """A missing required tool should appear in `warnings` (machine-readable)."""
    import json as _json

    runner = CliRunner()
    fake_env = {
        "ppanggolin": {"found": True, "version": "2.2", "path": "/x"},
        "gapseq": {"found": False, "version": None, "path": None},  # MISSING
        "diamond": {"found": True, "version": "2.1", "path": "/x"},
        "docker": {"found": True, "version": "24", "path": "/x"},
        "podman": {"found": False, "version": None, "path": None},
        "apptainer": {"found": False, "version": None, "path": None},
    }
    fake_data = {
        n: {"found": True, "size_bytes": 1, "rows": 1, "path": "/x"}
        for n in ("mnxr_to_rhea.tsv", "seed_to_rhea.tsv", "seed_to_bigg.tsv")
    }

    with (
        patch("parmmosahn.utils.check_environment", return_value=fake_env),
        patch("parmmosahn.utils.check_bundled_data", return_value=fake_data),
    ):
        result = runner.invoke(cli, ["doctor", "--json"])

    payload = _json.loads(result.output)
    assert payload["checks"]["all_required_tools_ok"] is False
    assert any("required external tools" in w for w in payload["warnings"])


def test_cli_doctor_warns_on_missing_data_file():
    """doctor surfaces a WARNING when a bundled TSV is unreachable."""
    runner = CliRunner()

    env_ok = {
        name: {"found": True, "version": "x", "path": "/x"}
        for name in ("ppanggolin", "gapseq", "diamond", "docker", "podman", "apptainer")
    }
    fake_data = {
        "mnxr_to_rhea.tsv": {
            "found": True,
            "size_bytes": 100,
            "rows": 10,
            "path": "/x",
        },
        "seed_to_rhea.tsv": {
            "found": True,
            "size_bytes": 100,
            "rows": 10,
            "path": "/x",
        },
        "seed_to_bigg.tsv": {
            "found": False,
            "size_bytes": None,
            "rows": None,
            "path": None,
        },
    }

    with (
        patch("parmmosahn.utils.check_environment", return_value=env_ok),
        patch("parmmosahn.utils.check_bundled_data", return_value=fake_data),
    ):
        result = runner.invoke(cli, ["doctor"])

    assert result.exit_code == 0
    assert "[MISSING] seed_to_bigg.tsv" in result.output
    assert "Annotation enrichment (Step 5b) will fail" in result.output


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


def test_cli_init_config_template_is_complete_and_parseable(tmp_path: Path):
    """Generated template must mention every config-loadable flag the project
    and pan commands honour, and must round-trip through YAML parsing."""
    pytest_importorskip_yaml()
    import yaml

    runner = CliRunner()
    out_file = tmp_path / "config.yml"
    runner.invoke(cli, ["init-config", "-o", str(out_file)])

    content = out_file.read_text()

    # Every key that `pan` and `project` look up via cfg.get must appear in
    # the template — either uncommented (defaults) or commented (examples).
    expected_keys = [
        "genomes",
        "fasta",
        "panmodel_rxns",
        "medium",
        "output",
        "label",
        "add_unique",
        "soft_core",
        "threads",
        "diamond_bits",
        "gapseq_bits",
        "biomass",
        "parallel_strains",
        "polish",
        "engine",
        "curation_db",
        "pause_for_curation",
    ]
    missing = [k for k in expected_keys if k not in content]
    assert not missing, f"Template is missing keys: {missing}"

    # The non-commented portion must parse as valid YAML on its own.
    uncommented = "\n".join(
        line for line in content.splitlines() if not line.lstrip().startswith("#")
    )
    parsed = yaml.safe_load(uncommented)
    # Defaults that the template sets explicitly:
    assert parsed["threads"] == 14
    assert parsed["diamond_bits"] == 150
    assert parsed["gapseq_bits"] == 150
    assert parsed["biomass"] == "pos"
    assert parsed["polish"] is False
    assert parsed["engine"] == "docker"


def test_cli_project_config_loads_downstream_keys(tmp_path: Path):
    """v0.3 init-config keys (parallel_strains, polish, curation_db,
    pause_for_curation) must reach derive_strain_models / run_downstream
    when the user supplies them via --config rather than CLI flags."""
    pytest_importorskip_yaml()

    genomes = tmp_path / "g"
    genomes.mkdir()
    fasta = tmp_path / "f.faa"
    fasta.touch()
    rxns = tmp_path / "r.tbl"
    rxns.touch()
    medium = tmp_path / "m.csv"
    medium.touch()
    curation_db = tmp_path / "curation.xlsx"
    curation_db.touch()
    out = tmp_path / "out"

    config_file = tmp_path / "config.yml"
    config_file.write_text(
        f"genomes: {genomes}\n"
        f"fasta: {fasta}\n"
        f"panmodel_rxns: {rxns}\n"
        f"medium: {medium}\n"
        f"curation_db: {curation_db}\n"
        f"parallel_strains: 3\n"
        f"polish: true\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    with (
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.pipeline.run_downstream") as mock_downstream,
    ):
        mock_derive.return_value = out / "03-strain_models" / "models"
        mock_downstream.return_value = {}
        result = runner.invoke(
            cli, ["project", "-o", str(out), "--config", str(config_file)]
        )

    assert result.exit_code == 0, result.output
    _, derive_kwargs = mock_derive.call_args
    assert derive_kwargs["parallel_strains"] == 3

    _, downstream_kwargs = mock_downstream.call_args
    # polish=True in config means no_polish=not polish=False reaches downstream
    assert downstream_kwargs["no_polish"] is False
    assert downstream_kwargs["curation_db"] == curation_db


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
        mock_derive.return_value = out / "03-strain_models" / "models"
        result = runner.invoke(
            cli,
            ["project", "-o", str(out), "--config", str(config_file)],
        )

    assert result.exit_code == 0, result.output
    mock_derive.assert_called_once()
    _, kwargs = mock_derive.call_args
    assert kwargs["threads"] == 6
    assert kwargs["diamond_bits"] == 111


def test_cli_run_config_loads_all_keys(tmp_path: Path):
    """'run --config' should pull every config-template key into the pipeline."""
    pytest_importorskip_yaml()

    genomes = tmp_path / "g"
    genomes.mkdir()
    medium = tmp_path / "m.csv"
    medium.touch()
    curation_db = tmp_path / "curation.xlsx"
    curation_db.touch()
    out = tmp_path / "out"

    config_file = tmp_path / "config.yml"
    config_file.write_text(
        f"genomes: {genomes}\n"
        f"output: {out}\n"
        f"label: cohort_test\n"
        f"medium: {medium}\n"
        f"threads: 6\n"
        f"diamond_bits: 200\n"
        f"gapseq_bits: 180\n"
        f"biomass: neg\n"
        f"parallel_strains: 4\n"
        f"polish: true\n"
        f"engine: podman\n"
        f"curation_db: {curation_db}\n"
        f"add_unique: true\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    with (
        patch("parmmosahn.pipeline.run_pan") as mock_run_pan,
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.pipeline.run_downstream") as mock_downstream,
        patch("parmmosahn.utils.validate_input_genomes"),
    ):
        # run_pan returns (step_times, rxn_table_path, softcore_fasta)
        mock_run_pan.return_value = (tmp_path / "rxns.tbl", tmp_path / "pan.faa")
        mock_derive.return_value = out / "03-strain_models" / "models"
        mock_downstream.return_value = {}

        result = runner.invoke(cli, ["run", "--config", str(config_file)])

    assert result.exit_code == 0, result.output

    # run_pan should see threads + add_unique from config
    _, pan_kwargs = mock_run_pan.call_args
    assert pan_kwargs["threads"] == 6
    assert pan_kwargs["add_unique"] is True

    # derive_strain_models should see parallel_strains + diamond_bits + biomass
    _, derive_kwargs = mock_derive.call_args
    assert derive_kwargs["parallel_strains"] == 4
    assert derive_kwargs["diamond_bits"] == 200
    assert derive_kwargs["biomass"] == "neg"

    # run_downstream should see polish, engine, curation_db
    _, downstream_kwargs = mock_downstream.call_args
    assert downstream_kwargs["no_polish"] is False  # polish=True → no_polish=False
    assert downstream_kwargs["engine"] == "podman"
    assert downstream_kwargs["curation_db"] == curation_db


def test_cli_run_cli_flag_overrides_config(tmp_path: Path):
    """An explicit --threads on the CLI must beat the same key in config."""
    pytest_importorskip_yaml()

    genomes = tmp_path / "g"
    genomes.mkdir()
    medium = tmp_path / "m.csv"
    medium.touch()

    config_file = tmp_path / "config.yml"
    config_file.write_text(
        f"genomes: {genomes}\n"
        f"output: {tmp_path / 'out'}\n"
        f"label: x\n"
        f"medium: {medium}\n"
        f"threads: 6\n",  # config says 6
        encoding="utf-8",
    )

    runner = CliRunner()
    with (
        patch("parmmosahn.pipeline.run_pan") as mock_run_pan,
        patch("parmmosahn.models.derive_strain_models") as mock_derive,
        patch("parmmosahn.pipeline.run_downstream"),
        patch("parmmosahn.utils.validate_input_genomes"),
    ):
        mock_run_pan.return_value = (tmp_path / "rxns.tbl", tmp_path / "pan.faa")
        mock_derive.return_value = tmp_path / "models"

        # CLI says 12 — must win.
        result = runner.invoke(cli, ["run", "--config", str(config_file), "-t", "12"])

    assert result.exit_code == 0, result.output
    _, pan_kwargs = mock_run_pan.call_args
    assert pan_kwargs["threads"] == 12


def test_cli_run_missing_required_args_with_no_config(tmp_path: Path):
    """`run` with neither CLI nor config should raise a clear UsageError."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run"])
    assert result.exit_code != 0
    assert "Missing required arguments" in result.output or "Usage" in result.output


def pytest_importorskip_yaml():
    """Skip test if pyyaml is not installed."""
    import pytest

    pytest.importorskip(
        "yaml", reason="pyyaml not installed (pip install parmmosahn[config])"
    )
