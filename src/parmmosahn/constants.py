"""Pipeline step labels and banner utilities for consistent log output."""

# Keyed by int for the seven numbered stages, plus string keys for the
# sub-steps (4a curation template, 5b annotation enrichment) so step_banner()
# is the single source of truth for every banner the pipeline emits.
STEPS: dict[int | str, str] = {
    1: "Building Pangenome",
    2: "Building Pan-model",
    3: "Deriving Strain Models",
    4: "Applying Automatic Curations",
    "4a": "Generating Curation Template",
    5: "Annotation Enrichment",
    "5b": "ModelPolisher",
    6: "Gathering Final Results",
    7: "Final Memote Evaluation",
}


def step_banner(step: int | str, detail: str = "") -> str:
    """Return a consistent step log banner string.

    Args:
        step: Pipeline step number or sub-step key (e.g. 1, "4a", "5b").
        detail: Optional parenthetical suffix (e.g. 'SKIPPED by user').

    Returns:
        Formatted banner string, e.g. '=== STEP 1: Building Pangenome ==='.
    """
    label = STEPS.get(step, f"Step {step}")
    suffix = f" ({detail})" if detail else ""
    return f"=== STEP {step}: {label}{suffix} ==="
