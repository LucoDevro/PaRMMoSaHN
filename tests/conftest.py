import shutil
from pathlib import Path

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Biological test fixtures (session-scoped, point to test_run/ artifacts)
# ---------------------------------------------------------------------------

_TEST_RUN = Path(__file__).resolve().parent.parent / "test_run"

EXPECTED_STRAINS = [
    "GCA_000227485.1",
    "GCF_000009045.1",
    "GCF_000146565.1",
    "GCF_000209795.2",
]


@pytest.fixture(scope="session")
def test_run_dir() -> Path:
    if not _TEST_RUN.is_dir():
        pytest.skip("test_run/ directory not present")
    return _TEST_RUN


@pytest.fixture(scope="session")
def precomputed_pan_fasta(test_run_dir: Path) -> Path:
    p = test_run_dir / "01-pangenome" / "test_run.faa"
    if not p.exists():
        pytest.skip("Pre-computed pan FASTA not found")
    return p


@pytest.fixture(scope="session")
def precomputed_pan_rxn_table(test_run_dir: Path) -> Path:
    p = test_run_dir / "02-panmodel" / "test_run-all-Reactions.tbl"
    if not p.exists():
        pytest.skip("Pre-computed pan reaction table not found")
    return p


@pytest.fixture(scope="session")
def strain_output_root(test_run_dir: Path) -> Path:
    """Return the root directory containing matches/, proteomes/, models/."""
    for candidate in ["03-strain_models", "output/03-strains"]:
        root = test_run_dir / candidate
        if (root / "models").is_dir():
            return root
    pytest.skip("No strain output directory found in test_run/")
    # Satisfy static type checkers: pytest.skip raises, but type checkers
    # may not recognise that as a NoReturn. Return a Path as a fallback.
    return test_run_dir


@pytest.fixture(scope="session")
def available_strains(strain_output_root: Path) -> list[str]:
    """Return strains that have completed model directories."""
    models_root = strain_output_root / "models"
    present = [s for s in EXPECTED_STRAINS if (models_root / s).is_dir()]
    if not present:
        pytest.skip("No completed strain models found")
    return present


@pytest.fixture(scope="session")
def precomputed_strain_dirs(
    strain_output_root: Path, available_strains: list[str]
) -> dict[str, Path]:
    models_root = strain_output_root / "models"
    return {s: models_root / s for s in available_strains}


@pytest.fixture(scope="session")
def loaded_cobra_models(precomputed_strain_dirs: dict[str, Path]):
    import cobra as cb

    models = {}
    for strain, d in precomputed_strain_dirs.items():
        xml = d / f"{strain}.xml"
        if xml.exists():
            models[strain] = cb.io.read_sbml_model(str(xml))
    if not models:
        pytest.skip("No XML models found to load")
    return models


@pytest.fixture(scope="session")
def strain_genomes_dir(test_run_dir: Path) -> Path:
    d = test_run_dir / "genomes"
    if not d.is_dir():
        pytest.skip("test_run/genomes/ not found")
    return d


@pytest.fixture(scope="session")
def medium_file(test_run_dir: Path) -> Path:
    p = test_run_dir / "M9_glc_aerobic.csv"
    if not p.exists():
        pytest.skip("Medium file not found")
    return p


@pytest.fixture(scope="session")
def require_pipeline_tools():
    """Skip tests if external pipeline tools are not installed."""
    from parmmosahn.utils import check_environment

    results = check_environment()
    missing = [t for t in ["diamond", "gapseq"] if not results[t]["found"]]
    if missing:
        pytest.skip(f"External tools not available: {', '.join(missing)}")


# ---------------------------------------------------------------------------
# Unit test fixtures
# ---------------------------------------------------------------------------


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
