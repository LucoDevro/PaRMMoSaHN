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


# ---------------------------------------------------------------------------
# 4. Pan-reactome Characterization
# ---------------------------------------------------------------------------


def panreactome(
    models_dir: Path,
    output_path: Path,
    plot: bool = False,
) -> Path:
    """Characterize the pan-reactome across all strain models.

    Produces a multi-file output:
    - ``{output_path}`` (TSV): per-reaction classification (core/accessory/unique)
    - ``{stem}_summary.tsv``: pan-reactome size statistics
    - ``{stem}_jaccard.tsv``: pairwise Jaccard similarity matrix
    - ``{stem}_accumulation.tsv``: reaction accumulation curve
    - ``{stem}_dendrogram.png``: hierarchical clustering (if *plot* is True)

    Parameters
    ----------
    models_dir:
        Directory containing SBML ``.xml`` files.
    output_path:
        Path for the main output TSV file.
    plot:
        If *True*, save a dendrogram PNG showing strain clustering by
        metabolic similarity.  Requires ``matplotlib`` (``pip install
        'parmmosahn[viz]'``).

    Returns
    -------
    Path to the main reaction classification TSV.
    """
    models = _load_models(models_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_strains = len(models)

    strain_names = [name for name, _ in models]

    # Build reaction presence sets
    rxn_sets: dict[str, set[str]] = {}
    all_rxns: set[str] = set()
    for name, model in models:
        rxns = {r.id for r in model.reactions}
        rxn_sets[name] = rxns
        all_rxns.update(rxns)

    rxn_list = sorted(all_rxns)

    # --- Classify reactions ---
    rxn_class: dict[str, dict[str, Any]] = {}
    for rxn in rxn_list:
        present_in = [s for s in strain_names if rxn in rxn_sets[s]]
        count = len(present_in)
        if count == n_strains:
            category = "core"
        elif count == 1:
            category = "unique"
        else:
            category = "accessory"
        rxn_class[rxn] = {
            "count": count,
            "fraction": count / n_strains,
            "category": category,
            "strains": ";".join(present_in),
        }

    # Write main classification TSV
    rows: list[str] = ["reaction\tcount\tfraction\tcategory\tstrains"]
    for rxn in rxn_list:
        r = rxn_class[rxn]
        rows.append(
            f"{rxn}\t{r['count']}\t{r['fraction']:.4f}\t{r['category']}\t{r['strains']}"
        )
    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    # --- Summary statistics ---
    core_count = sum(1 for r in rxn_class.values() if r["category"] == "core")
    accessory_count = sum(1 for r in rxn_class.values() if r["category"] == "accessory")
    unique_count = sum(1 for r in rxn_class.values() if r["category"] == "unique")

    model_sizes = [len(rxn_sets[s]) for s in strain_names]
    mean_size = sum(model_sizes) / len(model_sizes)
    std_size = (
        sum((x - mean_size) ** 2 for x in model_sizes) / len(model_sizes)
    ) ** 0.5

    summary_path = output_path.with_name(f"{output_path.stem}_summary.tsv")
    summary_rows = [
        "metric\tvalue",
        f"n_strains\t{n_strains}",
        f"pan_reactome_size\t{len(rxn_list)}",
        f"core_reactions\t{core_count}",
        f"accessory_reactions\t{accessory_count}",
        f"unique_reactions\t{unique_count}",
        f"core_fraction\t{core_count / len(rxn_list):.4f}",
        f"mean_model_size\t{mean_size:.1f}",
        f"std_model_size\t{std_size:.1f}",
        f"min_model_size\t{min(model_sizes)}",
        f"max_model_size\t{max(model_sizes)}",
    ]
    summary_path.write_text("\n".join(summary_rows) + "\n", encoding="utf-8")

    logger.info(
        f"Pan-reactome: {len(rxn_list)} total reactions "
        f"({core_count} core, {accessory_count} accessory, {unique_count} unique) "
        f"across {n_strains} strains"
    )

    # --- Jaccard similarity matrix ---
    jaccard_path = output_path.with_name(f"{output_path.stem}_jaccard.tsv")
    jac_rows: list[str] = ["strain\t" + "\t".join(strain_names)]
    for s1 in strain_names:
        vals: list[str] = []
        for s2 in strain_names:
            intersection = len(rxn_sets[s1] & rxn_sets[s2])
            union = len(rxn_sets[s1] | rxn_sets[s2])
            jac = intersection / union if union > 0 else 0.0
            vals.append(f"{jac:.4f}")
        jac_rows.append(f"{s1}\t" + "\t".join(vals))
    jaccard_path.write_text("\n".join(jac_rows) + "\n", encoding="utf-8")

    # --- Accumulation curve ---
    import random

    accum_path = output_path.with_name(f"{output_path.stem}_accumulation.tsv")
    # Average over multiple random orderings for a smoother curve
    n_permutations = min(100, max(10, n_strains * 5))
    accum_sums = [0.0] * n_strains

    rng = random.Random(42)  # deterministic seed
    for _ in range(n_permutations):
        order = list(range(n_strains))
        rng.shuffle(order)
        seen: set[str] = set()
        for pos, idx in enumerate(order):
            seen.update(rxn_sets[strain_names[idx]])
            accum_sums[pos] += len(seen)

    accum_rows = ["n_genomes\tmean_pan_reactions"]
    for i in range(n_strains):
        accum_rows.append(f"{i + 1}\t{accum_sums[i] / n_permutations:.1f}")
    accum_path.write_text("\n".join(accum_rows) + "\n", encoding="utf-8")

    logger.info(
        f"Pan-reactome analysis written to {output_path.parent}: "
        f"classification, summary, Jaccard matrix, accumulation curve"
    )

    if plot:
        _plot_dendrogram(jaccard_path, strain_names, rxn_sets)

    return output_path


def _plot_dendrogram(
    jaccard_tsv: Path,
    strain_names: list[str],
    rxn_sets: dict[str, set[str]],
) -> Path:
    """Render a dendrogram from Jaccard distances between strain metabolic models."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform
    except ImportError:
        logger.warning(
            "matplotlib/scipy not installed — skipping dendrogram. "
            "Install with: pip install 'parmmosahn[viz]' scipy"
        )
        return jaccard_tsv

    n = len(strain_names)
    dist_matrix = np.zeros((n, n))
    for i, s1 in enumerate(strain_names):
        for j, s2 in enumerate(strain_names):
            if i != j:
                intersection = len(rxn_sets[s1] & rxn_sets[s2])
                union = len(rxn_sets[s1] | rxn_sets[s2])
                dist_matrix[i, j] = 1.0 - (intersection / union if union > 0 else 0.0)

    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method="ward")

    png_path = jaccard_tsv.with_name(
        jaccard_tsv.stem.replace("_jaccard", "_dendrogram") + ".png"
    )
    fig, ax = plt.subplots(figsize=(max(6, n * 0.8), 5))
    dendrogram(Z, labels=strain_names, ax=ax, leaf_rotation=90)
    ax.set_ylabel("Jaccard distance (1 - similarity)")
    ax.set_title("Strain clustering by metabolic reaction content")
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)

    logger.info(f"Dendrogram saved to {png_path}")
    return png_path


# ---------------------------------------------------------------------------
# 5. Phenotype Validation
# ---------------------------------------------------------------------------


def validate_phenotypes(
    models_dir: Path,
    phenotype_path: Path,
    output_path: Path,
    medium_path: Path | None = None,
    threshold: float = 1e-3,
) -> Path:
    """Validate model predictions against experimental growth phenotypes.

    Parameters
    ----------
    models_dir:
        Directory containing SBML ``.xml`` files.
    phenotype_path:
        CSV/TSV with experimental growth data. Expected format:

        - **Matrix format**: first column = exchange reaction IDs (carbon
          sources), remaining columns = strain names, values = 1 (growth)
          or 0 (no growth).
        - **Long format**: columns ``strain``, ``carbon_source``, ``growth``
          (1/0 or True/False).

        Exchange reaction IDs should match the model (e.g., ``EX_glc__D_e0``).
    output_path:
        Path for the output TSV with predictions vs observations.
    medium_path:
        Optional base medium CSV.  For each carbon source test, this medium
        is loaded and the specific exchange reaction is toggled.
    threshold:
        Predicted growth rate below this counts as "no growth".

    Returns
    -------
    Path to the results TSV.
    """
    import polars as pl

    models = _load_models(models_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    strain_names = [name for name, _ in models]
    model_map = {name: model for name, model in models}

    # --- Parse phenotype data ---
    sep = "\t" if phenotype_path.suffix == ".tsv" else ","
    df = pl.read_csv(phenotype_path, separator=sep)
    cols = df.columns

    # Detect format: long (strain, carbon_source, growth) vs matrix
    long_cols = {"strain", "carbon_source", "growth"}
    if long_cols.issubset({c.lower() for c in cols}):
        # Normalize column names
        col_map = {c: c.lower() for c in cols}
        df = df.rename(col_map)
        phenotypes: list[tuple[str, str, bool]] = []
        for row in df.iter_rows(named=True):
            s = str(row["strain"])
            cs = str(row["carbon_source"])
            g = (
                bool(int(row["growth"]))
                if not isinstance(row["growth"], bool)
                else row["growth"]
            )
            phenotypes.append((s, cs, g))
    else:
        # Matrix format: first column = carbon_source, rest = strains
        phenotypes = []
        cs_col = cols[0]
        for row in df.iter_rows(named=True):
            cs = str(row[cs_col])
            for col in cols[1:]:
                val = row[col]
                if val is not None:
                    phenotypes.append((col, cs, bool(int(val))))

    # --- Run FBA predictions ---
    results: list[dict[str, Any]] = []
    tp = fp = tn = fn = 0

    for strain, carbon_source, observed_growth in phenotypes:
        if strain not in model_map:
            logger.debug(f"Skipping {strain}/{carbon_source}: strain not in models")
            continue

        model = model_map[strain]

        # Set up medium
        if medium_path:
            _apply_medium_csv(model, medium_path)

        medium = dict(model.medium)

        # Check if the exchange reaction exists
        if carbon_source not in model.reactions:
            results.append(
                {
                    "strain": strain,
                    "carbon_source": carbon_source,
                    "observed": observed_growth,
                    "predicted": False,
                    "growth_rate": 0.0,
                    "match": not observed_growth,
                    "note": "exchange_not_in_model",
                }
            )
            if observed_growth:
                fn += 1
            else:
                tn += 1
            continue

        # Enable the carbon source
        medium[carbon_source] = 1000.0
        model.medium = medium

        sol = model.optimize()
        growth_rate = sol.objective_value if sol.objective_value is not None else 0.0
        predicted_growth = growth_rate > threshold

        match = predicted_growth == observed_growth
        results.append(
            {
                "strain": strain,
                "carbon_source": carbon_source,
                "observed": observed_growth,
                "predicted": predicted_growth,
                "growth_rate": growth_rate,
                "match": match,
                "note": "",
            }
        )

        if predicted_growth and observed_growth:
            tp += 1
        elif predicted_growth and not observed_growth:
            fp += 1
        elif not predicted_growth and not observed_growth:
            tn += 1
        else:
            fn += 1

    # --- Write results ---
    rows: list[str] = [
        "strain\tcarbon_source\tobserved\tpredicted\tgrowth_rate\tmatch\tnote"
    ]
    for r in results:
        rows.append(
            f"{r['strain']}\t{r['carbon_source']}\t"
            f"{'1' if r['observed'] else '0'}\t{'1' if r['predicted'] else '0'}\t"
            f"{r['growth_rate']:.6f}\t{'1' if r['match'] else '0'}\t{r['note']}"
        )

    # Summary metrics
    total = tp + fp + tn + fn
    accuracy = (tp + tn) / total if total > 0 else 0.0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    # Matthews Correlation Coefficient
    mcc_denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = ((tp * tn) - (fp * fn)) / mcc_denom if mcc_denom > 0 else 0.0

    rows.append("")
    rows.append(f"# Summary: {total} tests, {tp} TP, {fp} FP, {tn} TN, {fn} FN")
    rows.append(f"# Accuracy: {accuracy:.4f}")
    rows.append(f"# Sensitivity: {sensitivity:.4f}")
    rows.append(f"# Specificity: {specificity:.4f}")
    rows.append(f"# MCC: {mcc:.4f}")

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    logger.info(
        f"Phenotype validation: {total} tests, accuracy={accuracy:.3f}, "
        f"sensitivity={sensitivity:.3f}, specificity={specificity:.3f}, MCC={mcc:.3f}"
    )

    return output_path


# ---------------------------------------------------------------------------
# 6. Gene Essentiality Prediction
# ---------------------------------------------------------------------------


def gene_essentiality(
    models_dir: Path,
    output_path: Path,
    medium_path: Path | None = None,
    threshold: float = 1e-3,
) -> Path:
    """Predict essential genes via single-gene deletion FBA.

    Parameters
    ----------
    models_dir:
        Directory containing SBML ``.xml`` files.
    output_path:
        Path for the output TSV file.
    medium_path:
        Optional CSV to override each model's medium.
    threshold:
        Growth rate below this fraction of wild-type marks a gene as essential.

    Returns
    -------
    Path to the written TSV.
    """
    from cobra.flux_analysis import single_gene_deletion

    models = _load_models(models_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[str] = ["strain\tgene\twt_growth\tko_growth\tessential"]

    for name, model in models:
        if medium_path:
            _apply_medium_csv(model, medium_path)

        wt_sol = model.optimize()
        wt_growth = (
            wt_sol.objective_value if wt_sol.objective_value is not None else 0.0
        )

        if wt_growth < threshold:
            logger.warning(
                f"{name}: wild-type growth is {wt_growth:.6f} (below threshold), "
                f"skipping essentiality analysis"
            )
            continue

        logger.info(
            f"  {name}: testing {len(model.genes)} genes (WT growth={wt_growth:.4f})..."
        )

        result = single_gene_deletion(model)

        for _, row in result.iterrows():
            gene_ids = row["ids"]
            ko_growth = row["growth"] if row["growth"] is not None else 0.0
            gene_id = next(iter(gene_ids)) if gene_ids else "unknown"
            essential = ko_growth < threshold * wt_growth
            rows.append(
                f"{name}\t{gene_id}\t{wt_growth:.6f}\t{ko_growth:.6f}\t"
                f"{'1' if essential else '0'}"
            )

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(f"Gene essentiality predictions written to {output_path}")
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
