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
    xml_models = list(model_dir.rglob("*.xml"))
    if not xml_models:
        logger.warning(f"No .xml models found in {model_dir} to gather.")
        return output_dir

    logger.info(f"Gathering and converting {len(xml_models)} models...")

    n = len(xml_models)
    failed_models = []
    for i, xml_path in enumerate(xml_models, 1):
        strain_name = xml_path.stem
        logger.debug(f"[{i}/{n}] Processing final model: {strain_name}")

        # 1. Copy the original SBML file
        dest_xml = output_dir / xml_path.name
        shutil.copy2(xml_path, dest_xml)

        try:
            # 2. Load the model using cobrapy
            model = cb.io.read_sbml_model(str(dest_xml))

            # 3. Convert and save to JSON
            dest_json = output_dir / f"{strain_name}.json"
            cb.io.save_json_model(model, str(dest_json))

            # 4. Convert and save to MATLAB (.mat)
            dest_mat = output_dir / f"{strain_name}.mat"
            cb.io.save_matlab_model(model, str(dest_mat))

        except Exception as e:
            logger.error(f"Failed to convert model {strain_name}: {e}")
            failed_models.append(strain_name)
            continue

    if failed_models:
        logger.warning(
            f"Conversion failed for {len(failed_models)} models: {', '.join(failed_models)}"
        )
    else:
        logger.info(f"All models successfully gathered and converted in {output_dir}")

    return output_dir
