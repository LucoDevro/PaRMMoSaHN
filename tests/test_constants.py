import pytest

from parmmosahn.constants import STEPS, step_banner


def test_steps_dict_covers_all_pipeline_steps():
    """STEPS should define labels for steps 1 through 7."""
    for i in range(1, 8):
        assert i in STEPS
        assert isinstance(STEPS[i], str)
        assert len(STEPS[i]) > 0


def test_step_banner_basic():
    """step_banner should produce the expected '=== STEP N: Label ===' format."""
    banner = step_banner(1)
    assert banner == "=== STEP 1: Building Pangenome ==="


def test_step_banner_with_detail():
    """step_banner with a detail arg should append it in parentheses."""
    banner = step_banner("5b", "docker")
    assert banner == "=== STEP 5b: ModelPolisher (docker) ==="


def test_step_banner_skip_detail():
    """step_banner(N, 'SKIPPED by user') should include the skip notice."""
    banner = step_banner("5b", "SKIPPED by user")
    assert "SKIPPED by user" in banner
    assert "STEP 5b" in banner


def test_step_banner_unknown_step():
    """step_banner should handle unknown step numbers gracefully."""
    banner = step_banner(99)
    assert "STEP 99" in banner
    assert "Step 99" in banner  # fallback label


def test_step_banner_no_detail_no_parens():
    """Without a detail, the banner should not contain parentheses."""
    banner = step_banner(3)
    assert "(" not in banner
    assert ")" not in banner
