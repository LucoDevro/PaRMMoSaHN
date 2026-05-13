import json
from pathlib import Path
from unittest.mock import patch

import pytest

from parmmosahn.pipeline import run_downstream, run_pan, write_pipeline_summary


@patch("parmmosahn.pipeline.evaluate_final_models")
@patch("parmmosahn.pipeline.gather_and_convert_models")
@patch("parmmosahn.pipeline.apply_curations")
@patch("parmmosahn.pipeline.build_pan_model")
@patch("parmmosahn.pipeline.build_pangenome")
def test_run_pan(
    mock_build_pangenome,
    mock_build_pan_model,
    mock_evaluate,
    mock_gather,
    mock_curations,
    tmp_path: Path,
):
    """run_pan should call build_pangenome then build_pan_model and return correct paths."""
    label = "clostridiales"
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    output_dir = tmp_path / "pan_out"

    dummy_fasta = tmp_path / "dummy.faa"
    dummy_fasta.touch()
    mock_build_pangenome.return_value = dummy_fasta

    panmodel_rxns, softcore_fasta = run_pan(genomes_dir, output_dir, label, threads=4)

    mock_build_pangenome.assert_called_once_with(
        genomes_dir,
        output_dir / "01-pangenome",
        label,
        threads=4,
        add_unique=False,
        soft_core_override=None,
    )
    mock_build_pan_model.assert_called_once_with(
        dummy_fasta, output_dir / "02-panmodel", threads=4
    )

    assert panmodel_rxns == output_dir / "02-panmodel" / f"{label}-all-Reactions.tbl"
    assert softcore_fasta == dummy_fasta


@patch("parmmosahn.pipeline.evaluate_final_models")
@patch("parmmosahn.pipeline.gather_and_convert_models")
@patch("parmmosahn.pipeline.apply_curations")
@patch("parmmosahn.pipeline.build_pan_model")
@patch("parmmosahn.pipeline.build_pangenome")
def test_run_pan_with_soft_core_override(
    mock_build_pangenome,
    mock_build_pan_model,
    mock_evaluate,
    mock_gather,
    mock_curations,
    tmp_path: Path,
):
    """run_pan should forward soft_core_override to build_pangenome."""
    dummy_fasta = tmp_path / "dummy.faa"
    dummy_fasta.touch()
    mock_build_pangenome.return_value = dummy_fasta

    run_pan(tmp_path / "g", tmp_path / "out", "test", threads=2, soft_core_override=0.5)

    _, kwargs = mock_build_pangenome.call_args
    assert kwargs["soft_core_override"] == 0.5


@patch("parmmosahn.pipeline.evaluate_final_models")
@patch("parmmosahn.pipeline.gather_and_convert_models")
@patch(
    "parmmosahn.pipeline.enrich_models_directory",
    side_effect=lambda src, dst, proteomes_dir=None: dst,
)
@patch("parmmosahn.pipeline.apply_curations")
@patch("shutil.copy")
def test_run_downstream_with_curation(
    mock_copy,
    mock_apply_curations,
    mock_enrich,
    mock_gather,
    mock_evaluate,
    tmp_path: Path,
):
    """run_downstream with a curation DB should call apply_curations on each draft model."""
    draft_dir = tmp_path / "drafts"
    draft_dir.mkdir()

    # Create a nested draft model (strain-specific subdirectory layout)
    strain_dir = draft_dir / "strain_A"
    strain_dir.mkdir()
    draft_xml = strain_dir / "strain_A.xml"
    draft_xml.touch()

    curation_db = tmp_path / "curation.xlsx"
    curation_db.touch()

    output_dir = tmp_path / "out"
    run_downstream(
        draft_models_dir=draft_dir,
        output_dir=output_dir,
        curation_db=curation_db,
        no_polish=True,
    )

    expected_curated_out = output_dir / "04-curated_models"
    # apply_curations now receives a `stats` accumulator dict so the pipeline
    # can emit a per-row applicability report after the loop.
    mock_apply_curations.assert_called_once()
    args, kwargs = mock_apply_curations.call_args
    assert args == (draft_xml, curation_db, expected_curated_out)
    assert "stats" in kwargs and isinstance(kwargs["stats"], dict)
    mock_gather.assert_called_once()
    mock_evaluate.assert_called_once()


@patch("parmmosahn.pipeline.evaluate_final_models")
@patch("parmmosahn.pipeline.gather_and_convert_models")
@patch(
    "parmmosahn.pipeline.enrich_models_directory",
    side_effect=lambda src, dst, proteomes_dir=None: dst,
)
@patch("parmmosahn.pipeline.apply_curations")
@patch("parmmosahn.pipeline.generate_curation_template")
def test_run_downstream_no_curation_generates_template(
    mock_gen_template,
    mock_apply_curations,
    mock_enrich,
    mock_gather,
    mock_evaluate,
    tmp_path: Path,
):
    """Without a curation DB, the pipeline should auto-generate a template and apply curations."""
    draft_dir = tmp_path / "drafts"
    draft_dir.mkdir()
    strain_dir = draft_dir / "strain_B"
    strain_dir.mkdir()
    draft_xml = strain_dir / "strain_B.xml"
    draft_xml.touch()

    # generate_curation_template returns the path to the generated Excel
    fake_db = tmp_path / "out" / "curation_template.xlsx"
    mock_gen_template.return_value = fake_db

    output_dir = tmp_path / "out"
    run_downstream(
        draft_models_dir=draft_dir,
        output_dir=output_dir,
        curation_db=None,
        no_polish=True,
    )

    mock_gen_template.assert_called_once_with(draft_dir, output_dir)
    expected_curated_out = output_dir / "04-curated_models"
    mock_apply_curations.assert_called_once()
    args, kwargs = mock_apply_curations.call_args
    assert args == (draft_xml, fake_db, expected_curated_out)
    assert "stats" in kwargs and isinstance(kwargs["stats"], dict)
    mock_gather.assert_called_once()


@patch("parmmosahn.pipeline.run_model_polisher")
@patch("parmmosahn.pipeline.evaluate_final_models")
@patch("parmmosahn.pipeline.gather_and_convert_models")
@patch(
    "parmmosahn.pipeline.enrich_models_directory",
    side_effect=lambda src, dst, proteomes_dir=None: dst,
)
@patch("shutil.copy")
@patch("shutil.which", return_value="/usr/bin/docker")
def test_run_downstream_polisher_called_when_engine_present(
    mock_which,
    mock_copy,
    mock_enrich,
    mock_gather,
    mock_evaluate,
    mock_polisher,
    tmp_path: Path,
):
    """ModelPolisher should be invoked when engine is found and no_polish=False."""
    draft_dir = tmp_path / "drafts"
    draft_dir.mkdir()
    polished_dir = tmp_path / "polished"
    mock_polisher.return_value = polished_dir

    run_downstream(
        draft_models_dir=draft_dir,
        output_dir=tmp_path / "out",
        no_polish=False,
        engine="docker",
    )

    mock_polisher.assert_called_once()


# ---------------------------------------------------------------------------
# write_pipeline_summary: A1 — run_environment + provenance enrichment
# ---------------------------------------------------------------------------


def test_write_pipeline_summary_includes_run_environment(tmp_path: Path):
    """Summary should embed host CPU/RAM/OS so a reviewer can reconstruct the run."""
    out = tmp_path / "out"
    out.mkdir()
    write_pipeline_summary(
        out, steps={"01-pangenome": {"status": "completed"}}, parameters={}
    )

    summary = json.loads((out / "pipeline_summary.json").read_text())
    assert "run_environment" in summary
    env = summary["run_environment"]
    # cpu_count and platform are present on every supported host
    assert "platform" in env
    assert "python_version" in env


def test_write_pipeline_summary_records_medium_sha(tmp_path: Path):
    """Provenance should pin the SHA-256 of the medium CSV used for the run."""
    out = tmp_path / "out"
    out.mkdir()
    medium = tmp_path / "medium.csv"
    medium.write_text("compounds,name,maxFlux\ncpd00001,water,1000\n")

    write_pipeline_summary(out, steps={}, parameters={"medium": str(medium)})
    summary = json.loads((out / "pipeline_summary.json").read_text())
    prov = summary["provenance"]
    assert prov["medium_path"] == str(medium)
    assert isinstance(prov["medium_sha256"], str)
    assert len(prov["medium_sha256"]) == 64  # hex SHA-256


def test_write_pipeline_summary_records_curation_db_sha(tmp_path: Path):
    """Provenance should also pin the curation DB SHA-256 when provided."""
    out = tmp_path / "out"
    out.mkdir()
    db = tmp_path / "curation.xlsx"
    db.write_bytes(b"\x50\x4b\x03\x04 fake xlsx header")

    write_pipeline_summary(out, steps={}, parameters={"curation_db": str(db)})
    summary = json.loads((out / "pipeline_summary.json").read_text())
    prov = summary["provenance"]
    assert prov["curation_db_path"] == str(db)
    assert isinstance(prov["curation_db_sha256"], str)


def test_write_pipeline_summary_pulls_pangenome_meta(tmp_path: Path):
    """Provenance should include the soft-core threshold actually used."""
    out = tmp_path / "out"
    pan_dir = out / "01-pangenome"
    pan_dir.mkdir(parents=True)
    (pan_dir / "pangenome_meta.json").write_text(
        json.dumps(
            {
                "n_genomes": 794,
                "soft_core_threshold": 0.00252,
                "soft_core_override_used": False,
                "partitions": "softcore",
                "cluster_identity": 0.8,
                "cluster_coverage": 0.6,
                "add_unique": False,
            }
        )
    )

    write_pipeline_summary(out, steps={}, parameters={})
    summary = json.loads((out / "pipeline_summary.json").read_text())
    assert summary["provenance"]["pangenome"]["n_genomes"] == 794
    assert summary["provenance"]["pangenome"]["soft_core_threshold"] == 0.00252
