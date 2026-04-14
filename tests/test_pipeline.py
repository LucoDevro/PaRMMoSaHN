from pathlib import Path
from unittest.mock import patch

import pytest

from parmmosahn.pipeline import run_downstream, run_pan


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
    mock_apply_curations.assert_called_once_with(
        draft_xml, curation_db, expected_curated_out
    )
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
    mock_apply_curations.assert_called_once_with(
        draft_xml, fake_db, expected_curated_out
    )
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
