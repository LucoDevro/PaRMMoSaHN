"""
Optional analysis extensions for PaRMMoSaHN models.

Provides comparative FBA summaries, auxotrophy screening, and reaction
presence/absence heatmaps.  All functions work on directories of SBML models
produced by the main pipeline and write plain TSV files.
"""

import logging
import random
import polars as pl
from pathlib import Path
from typing import Any, Iterator, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

import cobra as cb
from cobra.core import Model, Solution
from cobra.flux_analysis import single_gene_deletion

from parmmosahn.utils import safe_init_output_file

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_models(models_dir: Path) -> Iterator[tuple[str, Model]]:
    """Yield ``(name, model)`` for each SBML model in *models_dir*, one at a time.

    Streaming (rather than returning a list) keeps memory flat: a 1000-model
    cohort is never held in RAM all at once. Subdirectories are searched, and
    models that fail to parse are logged and skipped. Raises ``FileNotFoundError``
    if the directory contains no ``.xml`` files.

    Because this is a generator, callers must iterate exactly once and keep only
    the cheap per-model summary they need (reaction-ID sets, growth rates, ...)
    rather than retaining the ``Model`` objects.
    """
    xml_files = sorted(models_dir.rglob("*.xml"))
    if not xml_files:
        raise FileNotFoundError(f"No .xml models found in {models_dir}")

    for xml in xml_files:
        try:
            model = cb.io.read_sbml_model(str(xml))
        except Exception as e:
            logger.warning(f"Skipping {xml.name}: {e}")
            continue
        yield xml.stem, model


def _safe_optimize(model: Model, name: str) -> Optional[Solution]:
    """Optimise *model*, returning the solution or ``None`` on failure.

    A solver error on a single pathological model must not abort a whole-cohort
    sweep, so the exception is caught, logged, and reported as ``None`` for the
    caller to record as an error row.
    """
    try:
        return model.optimize()
    except Exception as exc:
        logger.warning(f"  {name}: FBA failed ({exc}); recording as error")
        return None


def _apply_medium_csv(model: Model, medium_path: Path) -> None:
    """Override a model's medium from a two-column CSV (exchange_id, flux)."""
    df = pl.read_csv(medium_path)
    cols = df.columns
    medium: dict[str, float] = {}
    dropped = 0
    for row in df.iter_rows(named=True):
        rxn_id = str(row[cols[0]])
        try:
            flux = float(row[cols[1]])
        except (TypeError, ValueError):
            logger.warning(
                f"Medium CSV: skipping {rxn_id} — invalid flux {row[cols[1]]!r}"
            )
            continue
        # Ensure the exchange reaction exists in this model
        if rxn_id in model.reactions:
            medium[rxn_id] = flux
        else:
            dropped += 1
    if dropped:
        logger.debug(
            f"Medium CSV: {dropped} component(s) absent from model "
            f"{getattr(model, 'id', '?')} — skipped"
        )
    model.medium = medium
    
    return None


_TRUTHY = {"1", "true", "t", "yes", "y", "+"}
_FALSY = {"0", "false", "f", "no", "n", "-"}


def _coerce_growth(val: Any) -> bool:
    """
    Coerce a phenotype 'growth' cell to bool. Accepts bool, int/float (0/1),
    or common string forms ('True'/'False', 'yes'/'no', '1'/'0').
    """
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)):
        return bool(int(val))
    if isinstance(val, str):
        v = val.strip().lower()
        if v in _TRUTHY:
            return True
        if v in _FALSY:
            return False
    raise ValueError(f"Cannot interpret growth value: {val!r}")


# ---------------------------------------------------------------------------
# 1. FBA Summary
# ---------------------------------------------------------------------------


def fba_summary(
    models_dir: Path,
    output_path: Path,
    medium_path: Path | None = None,
    force: bool = False,
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
    force:
        Boolean flag to overwrite any existing output.

    Returns
    -------
    Path to the written TSV.
    """
    safe_init_output_file(output_path, force = force)

    rows: list[str] = []
    header = "strain\tstatus\tgrowth_rate\tnum_active_exchanges\tactive_exchanges"
    rows.append(header)

    for name, model in _load_models(models_dir):
        if medium_path:
            _apply_medium_csv(model, medium_path)

        sol = _safe_optimize(model, name)
        if sol is None:
            rows.append(f"{name}\terror\tNA\t0\t")
            continue

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
    force: bool = False,
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
    force:
        Boolean flag to overwrite any existing output.

    Returns
    -------
    Path to the written TSV.
    """
    safe_init_output_file(output_path, force = force)

    # Stream models once. For each strain we record only floats (its baseline
    # growth and per-component knockout growth), never the Model object, so the
    # whole cohort is never held in memory.
    strain_names: list[str] = []
    baselines: dict[str, float] = {}
    all_components: set[str] = set()
    # knockout_growth[strain][component] = growth with that component removed.
    knockout_growth: dict[str, dict[str, float]] = {}

    for name, model in _load_models(models_dir):
        if medium_path:
            _apply_medium_csv(model, medium_path)
        strain_names.append(name)

        med = dict(model.medium)
        all_components.update(med)

        base_sol = _safe_optimize(model, name)
        baselines[name] = (
            base_sol.objective_value
            if base_sol is not None and base_sol.objective_value is not None
            else 0.0
        )

        ko: dict[str, float] = {}
        for comp in sorted(med):
            # Build a fresh medium each iteration (never mutate `med` in place):
            # a failed optimisation then cannot leak a knocked-out component into
            # the next test, and no explicit restore is needed.
            knocked = dict(med)
            knocked[comp] = 0.0
            model.medium = knocked
            sol = _safe_optimize(model, name)
            ko[comp] = (
                sol.objective_value
                if sol is not None and sol.objective_value is not None
                else 0.0
            )
        knockout_growth[name] = ko

    components = sorted(all_components)

    # Write TSV. A component absent from a strain's medium reports that strain's
    # baseline growth (knocking out something it never consumed changes nothing).
    rows: list[str] = []
    rows.append("component\t" + "\t".join(strain_names))
    for comp in components:
        vals = [
            f"{knockout_growth[s].get(comp, baselines[s]):.6f}" for s in strain_names
        ]
        rows.append(f"{comp}\t" + "\t".join(vals))

    # Summary row: baseline growth
    rows.append("_baseline\t" + "\t".join(f"{baselines[s]:.6f}" for s in strain_names))

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(f"Auxotrophy screen written to {output_path}")

    # Log predicted auxotrophies
    for comp in components:
        for name in strain_names:
            growth = knockout_growth[name].get(comp, baselines[name])
            if baselines[name] > threshold and growth < threshold:
                logger.info(f"  Predicted auxotrophy: {name} requires {comp}")

    return output_path


# ---------------------------------------------------------------------------
# 3. Reaction Presence/Absence Heatmap
# ---------------------------------------------------------------------------


def reaction_heatmap(
    models_dir: Path,
    output_path: Path,
    plot: Path | None = None,
    force: bool = False,
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
    force:
        Boolean flag to overwrite any existing output.

    Returns
    -------
    Path to the written TSV.
    """
    safe_init_output_file(output_path, force = force)

    # Stream models once, keeping only each strain's reaction-ID set (cheap)
    # instead of holding the whole cohort in memory.
    rxn_sets: dict[str, set[str]] = {}
    all_rxns: set[str] = set()
    for name, model in _load_models(models_dir):
        rxns = {r.id for r in model.reactions}
        rxn_sets[name] = rxns
        all_rxns.update(rxns)

    strain_names = list(rxn_sets.keys())
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
        _plot_heatmap(plot, strain_names, rxn_list, rxn_sets, force = force)

    return output_path


# ---------------------------------------------------------------------------
# 4. Pan-reactome Characterization
# ---------------------------------------------------------------------------

# TODO: to revisit
def analyse_panreactome(
    models_dir: Path,
    output_path: Path,
    plot: bool = False,
    core_threshold: float = 0.99,
    cloud_threshold: float = 0.05,
    force: bool = False,
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
    core_threshold:
        A reaction is *core* if it occurs in strictly more than this fraction
        of strains. Default ``0.99`` (present in >99% of strains), the
        reaction-level adaptation of the Tettelin core/accessory/unique
        framework used throughout the thesis. The strictly-greater comparison
        matches ``pmth-gene-partition`` on the gene side: for ``N=794``,
        ``0.99 * 794 = 786.06`` so core requires ``>= 787`` strains.
    cloud_threshold:
        A reaction is *unique* (cloud) if it occurs in strictly fewer than this
        fraction of strains. Default ``0.05`` (<5%). Everything between the two
        thresholds is *accessory* (5-99%).
    force:
        Boolean flag to overwrite any existing output.

    Returns
    -------
    Path to the main reaction classification TSV.
    """
    summary_path = output_path.with_name(f"{output_path.stem}_summary.tsv")
    jaccard_path = output_path.with_name(f"{output_path.stem}_jaccard.tsv")
    accum_path = output_path.with_name(f"{output_path.stem}_accumulation.tsv")

    for file in (output_path, summary_path, jaccard_path, accum_path):
        safe_init_output_file(file, force = force)

    # Stream models once, keeping only each strain's reaction-ID set so the
    # whole cohort is never resident in memory at 1000-model scale.
    rxn_sets: dict[str, set[str]] = {}
    all_rxns: set[str] = set()
    for name, model in _load_models(models_dir):
        rxns = {r.id for r in model.reactions}
        rxn_sets[name] = rxns
        all_rxns.update(rxns)

    strain_names = list(rxn_sets.keys())
    n_strains = len(strain_names)
    if n_strains == 0:
        raise FileNotFoundError(f"No usable models found in {models_dir}")

    rxn_list = sorted(all_rxns)

    # --- Classify reactions ---
    # Thesis convention (reaction-level Tettelin core/accessory/unique):
    # core present in >core_threshold of strains (strictly greater), unique
    # in <cloud_threshold, accessory in between. Strictly-greater on the core
    # side matches the gene-side pmth-gene-partition classifier.
    core_cut = core_threshold * n_strains
    cloud_cut = cloud_threshold * n_strains
    rxn_class: dict[str, dict[str, Any]] = {}
    for rxn in rxn_list:
        present_in = [s for s in strain_names if rxn in rxn_sets[s]]
        count = len(present_in)
        if count > core_cut:
            category = "core"
        elif count < cloud_cut:
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

    summary_rows = [
        "metric\tvalue",
        f"n_strains\t{n_strains}",
        f"pan_reactome_size\t{len(rxn_list)}",
        f"core_reactions\t{core_count}",
        f"accessory_reactions\t{accessory_count}",
        f"unique_reactions\t{unique_count}",
        f"core_fraction\t{core_count / len(rxn_list):.4f}",
        f"core_threshold\t{core_threshold}",
        f"cloud_threshold\t{cloud_threshold}",
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
        _plot_dendrogram(jaccard_path, strain_names)

    return output_path


def _plot_dendrogram(
    jaccard_tsv: Path,
    strain_names: list[str],
) -> Path:
    """Render a dendrogram from Jaccard distances between strain metabolic models.

    Reads the pairwise Jaccard matrix already written to *jaccard_tsv* rather
    than recomputing all O(N^2) set intersections a second time.
    """
    n = len(strain_names)
    # Reuse the similarities on disk: row order matches strain_names (panreactome
    # writes the matrix in the same order). distance = 1 - jaccard similarity.
    dist_matrix = np.zeros((n, n))
    data_lines = jaccard_tsv.read_text(encoding="utf-8").strip().splitlines()[1:]
    for i, line in enumerate(data_lines):
        sims = line.split("\t")[1:]  # column 0 is the strain name
        for j, sim in enumerate(sims):
            dist_matrix[i, j] = 1.0 - float(sim)

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
    force: bool = False
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
    force:
        Boolean flag to overwrite any existing output.

    Returns
    -------
    Path to the results TSV.
    """
    safe_init_output_file(output_path, force = force)

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
            phenotypes.append((s, cs, _coerce_growth(row["growth"])))
    else:
        # Matrix format: first column = carbon_source, rest = strains
        phenotypes = []
        cs_col = cols[0]
        for row in df.iter_rows(named=True):
            cs = str(row[cs_col])
            for col in cols[1:]:
                val = row[col]
                if val is not None:
                    phenotypes.append((col, cs, _coerce_growth(val)))

    # Group tests by strain so we can stream models once (never holding the
    # whole cohort in memory) and run all of a strain's tests while it is loaded.
    tests_by_strain: dict[str, list[tuple[str, bool]]] = {}
    for strain, carbon_source, observed_growth in phenotypes:
        tests_by_strain.setdefault(strain, []).append((carbon_source, observed_growth))

    # --- Run FBA predictions ---
    results: list[dict[str, Any]] = []
    tp = fp = tn = fn = 0
    seen: set[str] = set()

    for name, model in _load_models(models_dir):
        if name not in tests_by_strain:
            continue
        seen.add(name)

        # Capture the baseline medium once, before any per-test mutation, so a
        # carbon source enabled for one test never leaks into the next.
        if medium_path:
            _apply_medium_csv(model, medium_path)
        base_medium = dict(model.medium)

        for carbon_source, observed_growth in tests_by_strain[name]:
            # Check if the exchange reaction exists
            if carbon_source not in model.reactions:
                results.append(
                    {
                        "strain": name,
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

            # Reset to baseline, then enable just this carbon source for the test.
            medium = dict(base_medium)
            medium[carbon_source] = 1000.0
            model.medium = medium

            sol = _safe_optimize(model, name)
            growth_rate = (
                sol.objective_value
                if sol is not None and sol.objective_value is not None
                else 0.0
            )
            predicted_growth = growth_rate > threshold

            match = predicted_growth == observed_growth
            results.append(
                {
                    "strain": name,
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

    # Strains present in the phenotype table but absent from the model directory.
    for strain in tests_by_strain:
        if strain not in seen:
            logger.debug(f"Skipping {strain}: strain not in models")

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
    essential_fraction: float = 1e-2,
    min_wt_growth: float = 1e-3,
    force: bool = False,
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
    essential_fraction:
        Relative cutoff: a knockout whose growth falls below this fraction of
        the wild-type optimum is called essential. Defaults to ``1e-2`` (1 %),
        matching the auxotrophy cutoff and the Methods description.
    min_wt_growth:
        Absolute viability floor: strains whose wild-type optimum is below this
        flux are skipped (the model cannot grow, so essentiality is undefined).
        This is a separate concern from ``essential_fraction`` and is kept at
        the conservative ``1e-3`` numerical-tolerance default.
    force:
        Boolean flag to overwrite any existing output.

    Returns
    -------
    Path to the written TSV.
    """
    safe_init_output_file(output_path, force = force)

    rows: list[str] = ["strain\tgene\twt_growth\tko_growth\tessential"]

    for name, model in _load_models(models_dir):
        if medium_path:
            _apply_medium_csv(model, medium_path)

        wt_sol = _safe_optimize(model, name)
        wt_growth = (
            wt_sol.objective_value
            if wt_sol is not None and wt_sol.objective_value is not None
            else 0.0
        )

        if wt_growth < min_wt_growth:
            logger.warning(
                f"{name}: wild-type growth is {wt_growth:.6f} (below the "
                f"{min_wt_growth:g} viability floor), skipping essentiality analysis"
            )
            continue

        logger.info(
            f"  {name}: testing {len(model.genes)} genes (WT growth={wt_growth:.4f})..."
        )

        try:
            result = single_gene_deletion(model)
        except Exception as exc:
            logger.warning(
                f"  {name}: single-gene deletion failed ({exc}); skipping strain"
            )
            continue

        for _, row in result.iterrows():
            gene_ids = row["ids"]
            ko_growth = row["growth"] if row["growth"] is not None else 0.0
            gene_id = next(iter(gene_ids)) if gene_ids else "unknown"
            essential = ko_growth < essential_fraction * wt_growth
            rows.append(
                f"{name}\t{gene_id}\t{wt_growth:.6f}\t{ko_growth:.6f}\t"
                f"{'1' if essential else '0'}"
            )

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(f"Gene essentiality predictions written to {output_path}")
    
    return output_path


def _plot_heatmap(
    plot_path: Path,
    strain_names: list[str],
    rxn_list: list[str],
    rxn_sets: dict[str, set[str]],
    force: bool = False,
) -> Path:
    """Render a clustered heatmap PNG from the reaction matrix."""
    safe_init_output_file(plot_path)

    # Build numpy array
    data = np.array(
        [[1 if rxn in rxn_sets[s] else 0 for s in strain_names] for rxn in rxn_list],
        dtype=np.int8,
    )

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
    g.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(g.figure)

    logger.info(f"Heatmap saved to {plot_path}")
    
    return plot_path


# ---------------------------------------------------------------------------
# Dead-end metabolite + reaction report (B2 dry-run stub)
# ---------------------------------------------------------------------------


def _classify_dead_end(metabolite) -> tuple[int, int]:
    """Return (producer_count, consumer_count) for a metabolite.

    Accounts for reaction reversibility: a reaction with a negative lower
    bound is a potential producer-or-consumer in both stoichiometric
    directions. A metabolite with zero producers is an orphan source;
    zero consumers is an orphan sink; either condition makes every
    reaction touching it unable to carry steady-state flux.
    """
    producers = 0
    consumers = 0
    for rxn in metabolite.reactions:
        coef = rxn.metabolites[metabolite]
        lb, ub = rxn.lower_bound, rxn.upper_bound
        if ub > 0:
            if coef > 0:
                producers += 1
            elif coef < 0:
                consumers += 1
        if lb < 0:
            if coef > 0:
                consumers += 1
            elif coef < 0:
                producers += 1
                
    return producers, consumers


def find_dead_ends(
    models_dir: Path,
    output_path: Path,
    force: bool = False
) -> Path:
    """Report dead-end metabolites + unambiguously-blocked reactions per model.

    Output columns:

    ``strain, total_metabolites, dead_end_metabolites, dead_end_fraction,
    total_reactions, reactions_all_dead_end, prunable_fraction,
    growth_rate_baseline, fba_status``

    Returns the path to the written TSV.
    """
    safe_init_output_file(output_path, force = force)
    
    rows: list[str] = [
        "\t".join(
            [
                "strain",
                "total_metabolites",
                "dead_end_metabolites",
                "dead_end_fraction",
                "total_reactions",
                "reactions_all_dead_end",
                "prunable_fraction",
                "growth_rate_baseline",
                "fba_status",
            ]
        )
    ]

    # Suppress cobra's per-load INFO chatter while streaming models one at a time.
    cobra_logger = logging.getLogger("cobra")
    prev_level = cobra_logger.level
    cobra_logger.setLevel(logging.WARNING)

    try:
        for name, model in _load_models(models_dir):
            mets = list(model.metabolites)
            dead_end_ids: set[str] = set()
            for met in mets:
                prods, cons = _classify_dead_end(met)
                if prods == 0 or cons == 0:
                    dead_end_ids.add(met.id)

            n_total = len(mets)
            n_dead = len(dead_end_ids)
            dead_frac = (n_dead / n_total) if n_total else 0.0

            rxns = list(model.reactions)
            n_rxns = len(rxns)
            all_dead_rxns = sum(
                1
                for r in rxns
                if r.metabolites and all(m.id in dead_end_ids for m in r.metabolites)
            )
            prun_frac = (all_dead_rxns / n_rxns) if n_rxns else 0.0

            try:
                sol = model.optimize()
                growth = (
                    f"{sol.objective_value:.6f}"
                    if sol.objective_value is not None
                    else "NA"
                )
                status = sol.status
            except Exception as exc:
                growth = "NA"
                status = f"error: {exc}"

            rows.append(
                "\t".join(
                    [
                        name,
                        str(n_total),
                        str(n_dead),
                        f"{dead_frac:.4f}",
                        str(n_rxns),
                        str(all_dead_rxns),
                        f"{prun_frac:.4f}",
                        growth,
                        status,
                    ]
                )
            )
            logger.info(
                f"  {name}: {n_dead}/{n_total} dead-end mets "
                f"({dead_frac:.1%}), {all_dead_rxns}/{n_rxns} fully-dead-end rxns"
            )
    finally:
        cobra_logger.setLevel(prev_level)

    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(f"Dead-ends report written to {output_path}")
    
    return output_path
