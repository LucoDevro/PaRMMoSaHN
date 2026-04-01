"""
Optional analysis extensions for PaRMMoSaHN models.

Provides comparative FBA summaries, auxotrophy screening, and reaction
presence/absence heatmaps.  All functions work on directories of SBML models
produced by the main pipeline and write plain TSV files.  Plotting requires
the optional ``viz`` dependency group (matplotlib + seaborn).
"""

import logging
from pathlib import Path
from typing import Any

import cobra as cb
from cobra.core import Model

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_models(models_dir: Path) -> list[tuple[str, Model]]:
    """Load all SBML models from *models_dir*, searching subdirectories too."""
    xml_files = sorted(models_dir.rglob("*.xml"))
    if not xml_files:
        raise FileNotFoundError(f"No .xml models found in {models_dir}")

    loaded: list[tuple[str, Model]] = []
    for xml in xml_files:
        try:
            model = cb.io.read_sbml_model(str(xml))
            loaded.append((xml.stem, model))
        except Exception as e:
            logger.warning(f"Skipping {xml.name}: {e}")
    return loaded


def _apply_medium_csv(model: Model, medium_path: Path) -> None:
    """Override a model's medium from a two-column CSV (exchange_id, flux)."""
    import polars as pl

    df = pl.read_csv(medium_path)
    cols = df.columns
    medium: dict[str, float] = {}
    for row in df.iter_rows(named=True):
        rxn_id = str(row[cols[0]])
        flux = float(row[cols[1]])
        # Ensure the exchange reaction exists in this model
        if rxn_id in model.reactions:
            medium[rxn_id] = flux
    model.medium = medium


# ---------------------------------------------------------------------------
# 1. FBA Summary
# ---------------------------------------------------------------------------


def fba_summary(
    models_dir: Path,
    output_path: Path,
    medium_path: Path | None = None,
) -> Path:
    """Run FBA on every model and write a comparative TSV.

    Columns: strain, status, growth_rate, num_active_exchanges, active_exchanges.

    Parameters
    ----------
    models_dir:
        Directory containing SBML ``.xml`` files (may include subdirs).
    output_path:
        Path for the output TSV file.
    medium_path:
        Optional CSV to override each model's medium (cols: exchange_id, flux).
        When *None* the medium already stored in each model is used.

    Returns
    -------
    Path to the written TSV.
    """
    models = _load_models(models_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[str] = []
    header = "strain\tstatus\tgrowth_rate\tnum_active_exchanges\tactive_exchanges"
    rows.append(header)

    for name, model in models:
        if medium_path:
            _apply_medium_csv(model, medium_path)

        sol = model.optimize()
        status = sol.status
        growth = (
            f"{sol.objective_value:.6f}" if sol.objective_value is not None else "NA"
        )

        # Collect exchange reactions with non-zero flux
        active: list[str] = []
        if sol.status == "optimal":
            for rxn in model.exchanges:
                flux = sol.fluxes.get(rxn.id, 0.0)
                if abs(flux) > 1e-8:
                    active.append(f"{rxn.id}={flux:.4f}")
        active.sort()

        rows.append(f"{name}\t{status}\t{growth}\t{len(active)}\t{';'.join(active)}")
        logger.info(f"  {name}: {status}, growth={growth}")

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(f"FBA summary written to {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# 2. Auxotrophy Screening
# ---------------------------------------------------------------------------


def auxotrophy_screen(
    models_dir: Path,
    output_path: Path,
    medium_path: Path | None = None,
    threshold: float = 1e-3,
) -> Path:
    """Systematically knock out each medium component and check growth.

    Produces a TSV matrix: rows = medium components, columns = strains,
    values = growth rate with that component removed.  A value near zero
    indicates a predicted auxotrophy.

    Parameters
    ----------
    models_dir:
        Directory containing SBML ``.xml`` files.
    output_path:
        Path for the output TSV file.
    medium_path:
        Optional CSV to override each model's medium.
    threshold:
        Growth rate below this is reported as auxotrophic (for the summary
        column only; the raw rates are always written).

    Returns
    -------
    Path to the written TSV.
    """
    models = _load_models(models_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if medium_path:
        for _, model in models:
            _apply_medium_csv(model, medium_path)

    # Collect the union of all medium components across models
    all_components: set[str] = set()
    for _, model in models:
        all_components.update(model.medium.keys())
    components = sorted(all_components)

    strain_names = [name for name, _ in models]

    # Baseline growth
    baselines: dict[str, float] = {}
    for name, model in models:
        sol = model.optimize()
        baselines[name] = sol.objective_value if sol.objective_value else 0.0

    # Matrix: component -> {strain: growth_rate_without_component}
    matrix: dict[str, dict[str, float]] = {}

    for comp in components:
        matrix[comp] = {}
        for name, model in models:
            medium = dict(model.medium)
            if comp not in medium:
                # Component not in this strain's medium — use baseline
                matrix[comp][name] = baselines[name]
                continue

            # Knock out this component
            original_flux = medium[comp]
            medium[comp] = 0.0
            model.medium = medium

            sol = model.optimize()
            growth = sol.objective_value if sol.objective_value is not None else 0.0
            matrix[comp][name] = growth

            # Restore
            medium[comp] = original_flux
            model.medium = medium

    # Write TSV
    rows: list[str] = []
    rows.append("component\t" + "\t".join(strain_names))
    for comp in components:
        vals = [f"{matrix[comp].get(s, 0.0):.6f}" for s in strain_names]
        rows.append(f"{comp}\t" + "\t".join(vals))

    # Summary row: baseline growth
    rows.append("_baseline\t" + "\t".join(f"{baselines[s]:.6f}" for s in strain_names))

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(f"Auxotrophy screen written to {output_path}")

    # Log predicted auxotrophies
    for comp in components:
        for name in strain_names:
            if baselines[name] > threshold and matrix[comp][name] < threshold:
                logger.info(f"  Predicted auxotrophy: {name} requires {comp}")

    return output_path


# ---------------------------------------------------------------------------
# 3. Reaction Presence/Absence Heatmap
# ---------------------------------------------------------------------------


def reaction_heatmap(
    models_dir: Path,
    output_path: Path,
    plot: bool = False,
) -> Path:
    """Build a reaction presence/absence matrix across strains.

    Parameters
    ----------
    models_dir:
        Directory containing SBML ``.xml`` files.
    output_path:
        Path for the output TSV file (the matrix).
    plot:
        If *True*, also save a clustered heatmap as PNG next to the TSV.
        Requires ``matplotlib`` and ``seaborn`` (install with
        ``pip install 'parmmosahn[viz]'``).

    Returns
    -------
    Path to the written TSV.
    """
    models = _load_models(models_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    strain_names = [name for name, _ in models]

    # Collect reaction sets
    rxn_sets: dict[str, set[str]] = {}
    all_rxns: set[str] = set()
    for name, model in models:
        rxns = {r.id for r in model.reactions}
        rxn_sets[name] = rxns
        all_rxns.update(rxns)

    rxn_list = sorted(all_rxns)

    # Build matrix
    rows: list[str] = []
    rows.append("reaction\t" + "\t".join(strain_names))
    for rxn in rxn_list:
        vals = ["1" if rxn in rxn_sets[s] else "0" for s in strain_names]
        rows.append(f"{rxn}\t" + "\t".join(vals))

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(
        f"Reaction matrix ({len(rxn_list)} reactions x {len(strain_names)} strains) written to {output_path}"
    )

    if plot:
        _plot_heatmap(output_path, strain_names, rxn_list, rxn_sets)

    return output_path


def _plot_heatmap(
    tsv_path: Path,
    strain_names: list[str],
    rxn_list: list[str],
    rxn_sets: dict[str, set[str]],
) -> Path:
    """Render a clustered heatmap PNG from the reaction matrix."""
    try:
        import matplotlib

        matplotlib.use("Agg")  # Non-interactive backend
        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sns
    except ImportError:
        logger.warning(
            "matplotlib/seaborn not installed — skipping heatmap plot. "
            "Install with: pip install 'parmmosahn[viz]'"
        )
        return tsv_path

    # Build numpy array
    data = np.array(
        [[1 if rxn in rxn_sets[s] else 0 for s in strain_names] for rxn in rxn_list],
        dtype=np.int8,
    )

    png_path = tsv_path.with_suffix(".png")

    # For large matrices, use a clustermap with no individual row labels
    fig_height = max(6, min(len(rxn_list) * 0.01, 40))
    fig_width = max(4, len(strain_names) * 1.2)

    g = sns.clustermap(
        data,
        xticklabels=strain_names,
        yticklabels=False,
        cmap="YlOrRd",
        figsize=(fig_width, fig_height),
        dendrogram_ratio=(0.1, 0.05),
        cbar_kws={"label": "Present"},
        method="ward",
        metric="euclidean",
    )
    g.ax_heatmap.set_xlabel("Strain")
    g.ax_heatmap.set_ylabel(f"Reactions (n={len(rxn_list)})")
    g.figure.suptitle("Reaction Presence/Absence", y=1.02)
    g.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(g.figure)

    logger.info(f"Heatmap saved to {png_path}")
    return png_path
