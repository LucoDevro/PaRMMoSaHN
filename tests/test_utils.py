import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from parmmosahn.utils import ToolExecutionError, extract_proteome_from_gbff, run_command

# ---------------------------------------------------------------------------
# run_command — uses subprocess.Popen + threads (not subprocess.run)
# ---------------------------------------------------------------------------


def test_run_command_success(tmp_path: Path):
    """run_command should return captured stdout on success."""
    output = run_command(["echo", "hello"])
    assert output is not None
    assert "hello" in output


def test_run_command_failure():
    """run_command should raise an exception when the subprocess exits non-zero."""
    with pytest.raises(Exception, match="Tool failed with exit code"):
        run_command(["false"])


def test_run_command_with_file_output(tmp_path: Path):
    """run_command should write stdout/stderr to files and return None."""
    stdout_file = tmp_path / "out.log"
    stderr_file = tmp_path / "err.log"

    result = run_command(
        ["echo", "file_output_test"],
        stdout_file=stdout_file,
        stderr_file=stderr_file,
    )

    assert result is None
    assert stdout_file.exists()
    assert "file_output_test" in stdout_file.read_text()


def test_run_command_stderr_captured_on_failure():
    """Stderr content should appear in the raised exception message."""
    with pytest.raises(Exception):
        # 'ls' on a non-existent path writes to stderr and exits non-zero
        run_command(["ls", "/path/that/does/not/exist/12345abcde"])


# ---------------------------------------------------------------------------
# extract_proteome_from_gbff
# ---------------------------------------------------------------------------


def test_extract_proteome_from_gbff(tmp_path: Path, dummy_gbff: Path):
    """GBFF files should be correctly converted to FAA files."""
    out_faa = tmp_path / "test.faa"

    result_path = extract_proteome_from_gbff(dummy_gbff, out_faa)

    assert result_path == out_faa
    assert out_faa.exists()
    content = out_faa.read_text()
    assert content.startswith(">")
    assert len(content) > 10


def test_extract_proteome_from_gbff_skips_existing(tmp_path: Path, dummy_gbff: Path):
    """If the output FAA already exists, extraction should be skipped."""
    out_faa = tmp_path / "existing.faa"
    out_faa.write_text(">already_here\nACGT\n")

    result_path = extract_proteome_from_gbff(dummy_gbff, out_faa)

    assert result_path == out_faa
    # File content should be unchanged (no re-extraction)
    assert out_faa.read_text() == ">already_here\nACGT\n"


# ---------------------------------------------------------------------------
# check_environment
# ---------------------------------------------------------------------------


def test_check_environment_found_tools():
    """Tools found in PATH should have found=True and a version string."""
    from parmmosahn.utils import check_environment

    fake_proc = MagicMock()
    fake_proc.stdout = "sometool 1.2.3\n"
    fake_proc.stderr = ""

    with (
        patch("parmmosahn.utils.shutil.which", return_value="/usr/bin/sometool"),
        patch("parmmosahn.utils.subprocess.run", return_value=fake_proc),
    ):
        results = check_environment()

    for name, info in results.items():
        assert info["found"] is True
        assert info["version"] is not None
        assert info["path"] == "/usr/bin/sometool"


def test_check_environment_missing_tools():
    """Tools absent from PATH should have found=False."""
    from parmmosahn.utils import check_environment

    with patch("parmmosahn.utils.shutil.which", return_value=None):
        results = check_environment()

    for name, info in results.items():
        assert info["found"] is False
        assert info["version"] is None
        assert info["path"] is None


def test_check_environment_version_timeout():
    """Timeout during version probing should not crash; version should be 'unknown'."""
    from parmmosahn.utils import check_environment

    with (
        patch("parmmosahn.utils.shutil.which", return_value="/usr/bin/tool"),
        patch(
            "parmmosahn.utils.subprocess.run",
            side_effect=subprocess.TimeoutExpired(cmd="tool", timeout=30),
        ),
    ):
        results = check_environment()

    for info in results.values():
        assert info["found"] is True
        assert info["version"] == "unknown"


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------


def test_load_config_basic(tmp_path: Path):
    """load_config should parse a YAML file and return a dict."""
    pytest.importorskip("yaml", reason="pyyaml not installed")
    from parmmosahn.utils import load_config

    cfg_file = tmp_path / "config.yml"
    cfg_file.write_text("threads: 8\nbiomass: neg\n", encoding="utf-8")

    cfg = load_config(cfg_file)
    assert cfg["threads"] == 8
    assert cfg["biomass"] == "neg"


def test_load_config_empty_file(tmp_path: Path):
    """An empty YAML file should return an empty dict."""
    pytest.importorskip("yaml", reason="pyyaml not installed")
    from parmmosahn.utils import load_config

    cfg_file = tmp_path / "empty.yml"
    cfg_file.write_text("", encoding="utf-8")

    assert load_config(cfg_file) == {}


def test_load_config_missing_pyyaml(tmp_path: Path):
    """load_config should raise ImportError with an install hint when pyyaml is absent."""
    from parmmosahn.utils import load_config

    cfg_file = tmp_path / "config.yml"
    cfg_file.write_text("threads: 4\n", encoding="utf-8")

    with patch.dict("sys.modules", {"yaml": None}):
        with pytest.raises(ImportError, match="parmmosahn\\[config\\]"):
            load_config(cfg_file)


# ---------------------------------------------------------------------------
# _extract_version
# ---------------------------------------------------------------------------


def test_extract_version_picks_semver_line():
    """_extract_version should prefer a line containing a version number."""
    from parmmosahn.utils import _extract_version

    output = "some tool\nsome tool version 1.4.2\nUsage: ..."
    assert "1.4.2" in _extract_version(output)


def test_extract_version_falls_back_to_first_line():
    """_extract_version should return the first non-empty line when no version found."""
    from parmmosahn.utils import _extract_version

    output = "No version here\nSecond line"
    assert _extract_version(output) == "No version here"


def test_extract_version_empty_output():
    """_extract_version should return 'unknown' for empty output."""
    from parmmosahn.utils import _extract_version

    assert _extract_version("") == "unknown"
    assert _extract_version("   \n  \n") == "unknown"


# ---------------------------------------------------------------------------
# sha256_file + get_system_resources (A1 + A3 helpers)
# ---------------------------------------------------------------------------


def test_sha256_file_returns_hex_digest(tmp_path: Path):
    """sha256_file should compute a stable hex digest for an existing file."""
    from parmmosahn.utils import sha256_file

    target = tmp_path / "medium.csv"
    target.write_text("compounds,name,maxFlux\ncpd00001,water,1000\n")

    digest = sha256_file(target)
    assert isinstance(digest, str)
    assert len(digest) == 64
    # Stable across calls
    assert sha256_file(target) == digest


def test_sha256_file_missing_returns_none(tmp_path: Path):
    """sha256_file should return None for non-existent paths instead of raising."""
    from parmmosahn.utils import sha256_file

    assert sha256_file(tmp_path / "does_not_exist") is None


def test_get_system_resources_basic_fields():
    """get_system_resources should always return platform + python_version."""
    from parmmosahn.utils import get_system_resources

    info = get_system_resources()
    assert "platform" in info
    assert "python_version" in info
    # cpu_count is best-effort but always present on supported hosts
    if "cpu_count" in info:
        assert info["cpu_count"] >= 1


def test_check_bundled_data_finds_required_tsvs():
    """All three runtime TSV lookups must be reachable from an editable install."""
    from parmmosahn.utils import check_bundled_data

    result = check_bundled_data()
    for name in ("mnxr_to_rhea.tsv", "seed_to_rhea.tsv", "seed_to_bigg.tsv"):
        assert name in result, f"{name} not probed"
        info = result[name]
        assert info[
            "found"
        ], f"{name} missing — check pyproject.toml package-data and MANIFEST.in"
        assert info["rows"] is not None and info["rows"] > 0
        assert info["size_bytes"] is not None and info["size_bytes"] > 0


def test_get_system_resources_includes_ram_on_linux():
    """On Linux, /proc/meminfo is the source of total/available RAM."""
    from parmmosahn.utils import get_system_resources

    if not Path("/proc/meminfo").exists():
        pytest.skip("non-Linux host")

    info = get_system_resources()
    assert "total_ram_gb" in info
    assert info["total_ram_gb"] > 0
    assert "available_ram_gb" in info
