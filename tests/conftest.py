import shutil
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def test_data_dir() -> Path:
    """Provides the absolute path to the tests/data/ directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def dummy_sbml(test_data_dir: Path, tmp_path: Path) -> Path:
    """Provides a fresh, disposable copy of the tiny SBML model."""
    source = test_data_dir / "tiny_model.xml"
    dest = tmp_path / "tiny_model.xml"
    shutil.copy(source, dest)
    return dest


@pytest.fixture
def dummy_gbff(test_data_dir: Path, tmp_path: Path) -> Path:
    """Provides a fresh, disposable copy of the tiny GenBank genome."""
    source = test_data_dir / "tiny_genome.gbff"
    dest = tmp_path / "tiny_genome.gbff"
    shutil.copy(source, dest)
    return dest


@pytest.fixture
def dummy_curation_db(tmp_path: Path) -> Path:
    """
    Dynamically generates a mock Excel curation database with the required sheets
    so we don't have to commit a binary .xlsx file to git.
    """
    db_path = tmp_path / "curation_db.xlsx"

    # 1. Mock the 'corrected_balances' sheet (Removed the M_ prefixes)
    df_balances = pd.DataFrame(
        {"reaction_id": ["RXN_1", "RXN_MISSING"], "new_formula": ["A_c -> 2.0 B_c", ""]}
    )

    # 2. Mock the 'duplicate_reactions' sheet
    df_duplicates = pd.DataFrame(
        {
            "Pair": ["'RXN_DUPE_1', 'RXN_DUPE_2'"],
            "Decision": [
                1
            ],  # Decision '1' means remove the second reaction (RXN_DUPE_2)
        }
    )

    # Write them to the Excel file
    with pd.ExcelWriter(db_path) as writer:
        df_balances.to_excel(writer, sheet_name="corrected_balances", index=False)
        df_duplicates.to_excel(writer, sheet_name="duplicate_reactions", index=False)

    return db_path
