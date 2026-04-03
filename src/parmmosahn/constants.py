"""Pipeline step labels and banner utilities for consistent log output."""

STEPS: dict[int, str] = {
    1: "Building Pangenome",
    2: "Building Pan-model",
    3: "Deriving Strain Models",
    4: "Applying Automatic Curations",
    5: "ModelPolisher",
    6: "Gathering Final Results",
    7: "Final Memote Evaluation",
}


def step_banner(step: int, detail: str = "") -> str:
    """Return a consistent step log banner string.

    Args:
        step: Pipeline step number.
        detail: Optional parenthetical suffix (e.g. 'SKIPPED by user').

    Returns:
        Formatted banner string, e.g. '=== STEP 1: Building Pangenome ==='.
    """
    label = STEPS.get(step, f"Step {step}")
    suffix = f" ({detail})" if detail else ""
    return f"=== STEP {step}: {label}{suffix} ==="
