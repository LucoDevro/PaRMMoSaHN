import logging
import shutil
from pathlib import Path

import cobra as cb

logger = logging.getLogger(__name__)


def gather_and_convert_models(model_dir: str | Path, output_dir: str | Path) -> Path:
    """
    Copies final SBML models into a centralized output directory and
    converts them into JSON and MATLAB formats for broad compatibility.

    Args:
        model_dir: Directory containing the final annotated .xml models.
        output_dir: Destination directory for the final outputs.

    Returns:
        Path to the final output directory.
    """
    model_dir = Path(model_dir).resolve()
    output_dir = Path(output_dir).resolve()

    if not model_dir.exists():
        raise FileNotFoundError(f"Source model directory not found: {model_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all XML models
    all_xml = sorted(model_dir.rglob("*.xml"))
    if not all_xml:
        logger.warning(f"No .xml models found in {model_dir} to gather.")
        return output_dir

    # Collapse duplicate stems (e.g. gapseq trees keep {strain}.xml plus copies
    # in per-strain sub-dirs). Outputs are keyed by stem, so without this the
    # same strain would convert twice and the later path would silently
    # overwrite the earlier. First path wins.
    xml_models: list[Path] = []
    seen_stems: set[str] = set()
    for xml in all_xml:
        if xml.stem in seen_stems:
            continue
        seen_stems.add(xml.stem)
        xml_models.append(xml)
    if len(xml_models) < len(all_xml):
        logger.info(
            f"{len(all_xml)} .xml files collapse to {len(xml_models)} unique "
            f"strains (same-stem copies in sub-directories ignored)."
        )

    logger.info(f"Gathering and converting {len(xml_models)} models...")

    n = len(xml_models)
    failed_models = []
    for i, xml_path in enumerate(xml_models, 1):
        strain_name = xml_path.stem
        logger.debug(f"[{i}/{n}] Processing final model: {strain_name}")

        dest_xml = output_dir / xml_path.name
        dest_json = output_dir / f"{strain_name}.json"
        dest_mat = output_dir / f"{strain_name}.mat"

        try:
            # Convert from the validated SOURCE (not a freshly-made copy) so a
            # copy fault is never misattributed to conversion. Read first, so an
            # unreadable model leaves no stray .xml copy behind.
            model = cb.io.read_sbml_model(str(xml_path))
            shutil.copy2(xml_path, dest_xml)
            cb.io.save_json_model(model, str(dest_json))
            cb.io.save_matlab_model(model, str(dest_mat))

        except Exception as e:
            logger.error(f"Failed to convert model {strain_name}: {e}")
            failed_models.append(strain_name)
            # Don't leave a partial JSON/MAT pair for downstream consumers.
            for partial in (dest_json, dest_mat):
                partial.unlink(missing_ok=True)
            continue

    if failed_models:
        logger.warning(
            f"Conversion failed for {len(failed_models)} models: {', '.join(failed_models)}"
        )
    else:
        logger.info(f"All models successfully gathered and converted in {output_dir}")

    return output_dir
