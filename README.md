# PaRMMoSaHN

**Pa**ngenome **R**eference-based **M**etabolic **Mo**delling by **Sa**ving **H**omolog **N**etworks

[![CI](https://github.com/robbedewin/PaRMMoSaHN/workflows/CI/badge.svg)](https://github.com/robbedewin/PaRMMoSaHN/actions)
[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**PaRMMoSaHN** is a Python-based pipeline that bridges partitioned pangenome graphs with high-fidelity metabolic reconstruction.

Instead of building genome-scale metabolic models (GEMs) for hundreds of strains from scratch -- which is computationally expensive and produces inconsistently gap-filled models -- PaRMMoSaHN builds a single high-quality **pan-model** from the pangenome and rapidly derives strain-specific models by mapping individual strains to this reference via homology networks.

### Why this approach?

Existing tools either reconstruct each strain independently (gapseq, CarveMe) or require a manually curated species-specific reference (Bactabolize, KpSC pan). PaRMMoSaHN is the first **species-agnostic, automated pipeline** that leverages PPanGGOLiN pangenome partitioning to build a shared metabolic reference, then projects it onto strains using DIAMOND sequence homology. This produces models that are:

- **Consistent** -- all strains are derived from the same metabolic reference, eliminating artificial variation from independent gap-filling
- **Fast** -- strain projection takes minutes instead of hours per genome
- **Scalable** -- `.done` checkpointing and isolated evaluation processes handle hundreds of strains

## Pipeline

```
   Genomes (.gbff)
        |
        v
  [Step 1] PPanGGOLiN pangenome construction
        |   -> soft-core protein families FASTA
        v
  [Step 2] gapseq metabolic pathway prediction on pangenome
        |   -> pan-model reaction table
        v
  [Step 3] Per-strain model derivation
        |   DIAMOND blastp -> filter pan-model -> gapseq draft + gapfill
        |   -> one SBML model per strain
        v
  [Step 4] Automated curation (memote evaluation + duplicate/imbalance fixes)
        |   Optional: pause for manual curation via Excel spreadsheet
        v
  [Step 5] ModelPolisher (optional, requires Docker/Podman/Apptainer)
        |
  [Step 5b] Annotation enrichment (NCBI protein, BRENDA, RHEA)
        |
        v
  [Step 6] Gather & convert (SBML, JSON, MATLAB)
        |
        v
  [Step 7] Final memote FAIR-compliance evaluation
```

### Output structure

```
output/
├── 01-pangenome/                       # PPanGGOLiN pangenome + soft-core FASTA
│   └── pangenome_meta.json             # n_genomes, soft-core threshold, cluster params
├── 02-panmodel/                        # gapseq reaction/pathway tables
├── 03-strain_models/                   # Per-strain DIAMOND matches, proteomes, draft models
├── 04-curated_models/                  # Models after automated curation
│   └── curation_application_report.tsv # Per-row hit counts across strains (v0.2.1)
├── 05-polished_models/                 # ModelPolisher output (if enabled)
├── 05b-annotated_models/               # SBML enriched with NCBI/BRENDA/RHEA annotations
├── 06-final_models/                    # Final models in XML, JSON, and MATLAB formats
├── 07-memote_reports/                  # draft/ and final/ memote HTML+JSON reports
├── curation_template.xlsx              # Memote-derived spreadsheet for manual curation
├── pipeline_summary.json               # Run metadata, parameters, provenance, model statistics
├── run.log                             # Full INFO-level run log
└── errors.log                          # WARNING+ messages (only if errors occurred)
```

The `pipeline_summary.json` file includes a `provenance` block (SHA-256 of
the medium CSV and curation database, soft-core threshold actually used,
external tool versions) and a `run_environment` block (host CPU/RAM/OS), so
a reviewer can verify a run is reproducible without re-running anything.

## Installation

PaRMMoSaHN orchestrates several bioinformatics tools that cannot be installed via pip alone. Use Conda/Mamba to set up the environment.

### 1. Create the Conda environment

```bash
# mamba is recommended for faster dependency resolution
mamba env create -f environment.yml
conda activate parmmosahn_env
```

### 2. Install PaRMMoSaHN

```bash
pip install -e .
```

Or directly from GitHub:

```bash
pip install "git+https://github.com/robbedewin/PaRMMoSaHN.git"
```

### 3. Verify the installation

```bash
parmmosahn doctor
```

This checks that all required external tools (PPanGGOLiN, gapseq, DIAMOND)
and optional container engines (Docker, Podman, Apptainer) are available,
and reports host CPU/RAM with a recommended memote-worker count for the
`evaluate` step (helpful on memory-constrained hosts such as default WSL2
installations, where the default worker count can trigger
`BrokenProcessPool` errors).

## Quick Start

### Workflow A: Full automated pipeline

```bash
parmmosahn run \
  -g /path/to/genomes/ \
  -o ./results/ \
  -l clostridiales \
  -m medium.csv \
  -t 14 --parallel-strains 2
```

**Required inputs:**
| Option | Description |
|--------|-------------|
| `-g, --genomes` | Directory containing annotated genomes in GenBank format (`.gbff`) |
| `-o, --output` | Output base directory |
| `-l, --label` | Label for the pangenome (used in filenames) |
| `-m, --medium` | Growth medium CSV for gap-filling (gapseq format) |

**Optional parameters:**
| Option | Default | Description |
|--------|---------|-------------|
| `-t, --threads` | 75% of CPUs | Total CPU threads |
| `--parallel-strains` | auto (`threads // 4`) | Strains processed in parallel in Step 3; each worker gets `threads / N` CPU threads (peaks ~1.5 GB RAM each) |
| `--diamond-bits` | 150 | DIAMOND blastp bitscore threshold |
| `--gapseq-bits` | 150 | gapseq pathway search bitscore threshold |
| `--biomass` | `pos` | Biomass reaction type (`pos` or `neg`, matches gapseq's gram-stain templates) |
| `--add-unique` | off | Include unique/cloud genes (singletons) in the pan-model |
| `--soft-core` | 2/N | Override the soft-core frequency threshold (fraction in 0–1) |
| `--polish` | off | Enable ModelPolisher (requires Docker/Podman/Apptainer; see [ModelPolisher](#modelpolisher) below) |
| `-e, --engine` | `docker` | Container engine when `--polish` is enabled |
| `-c, --curation-db` | none | Path to an existing curation spreadsheet (skips auto-template generation) |
| `--pause-for-curation` | off | Pause after Step 3 for manual curation; resume with `parmmosahn project --resume` |

### Workflow B: Human-in-the-loop curation

For maximum model quality, pause the pipeline after draft model evaluation, manually review the curation spreadsheet, then resume:

```bash
# Step 1: Run pipeline and pause for manual curation
parmmosahn run \
  -g ./genomes/ -o ./results/ -l my_species -m medium.csv \
  --pause-for-curation

# -> Edit results/curation_template.xlsx in Excel
#    Fill the 'duplicate_reactions' and 'curated_imbalances' sheets

# Step 2: Resume pipeline with your curation decisions
parmmosahn project --resume -o ./results/
```

The curation template has three sheets:
- **duplicate_reactions** -- pairs of duplicate reactions with decision options (keep 1, keep 2, drop both, keep both)
- **curated_imbalances** -- mass/charge-imbalanced reactions with a column for corrected formulas
- **ignored_imbalances** -- reactions to leave intact despite imbalance (with justification)

### Workflow C: Rapid projection of new isolates

If you already have a pan-model and sequenced new strains, bypass the pangenome construction:

```bash
parmmosahn project \
  -g ./new_strains/ \
  -f ./results/01-pangenome/my_species.faa \
  -r ./results/02-panmodel/my_species-all-Reactions.tbl \
  -m medium.csv \
  -o ./projection/
```

## Analysis

Once you have strain-specific models, PaRMMoSaHN provides built-in analysis commands under `parmmosahn analyze` to explore metabolic diversity, validate predictions, and compare strains. These analyses are what make the models scientifically useful -- raw SBML files only become insights when you interrogate them.

### Pan-reactome characterization

The pan-reactome is the union of all metabolic reactions across your strains, analogous to the pangenome but at the metabolic level. Characterizing it reveals which metabolic capabilities are universally conserved (core), which are shared by subsets of strains (accessory), and which are strain-specific (unique). This is the central scientific output of a pangenome-scale metabolic study.

```bash
parmmosahn analyze panreactome \
  -M ./results/06-final_models/ \
  -o ./analysis/panreactome.tsv \
  --plot
```

**Outputs:**
- `panreactome.tsv` -- per-reaction classification (core/accessory/unique) with strain presence
- `panreactome_summary.tsv` -- pan-reactome size, core/accessory/unique counts, model size statistics
- `panreactome_jaccard.tsv` -- pairwise Jaccard similarity matrix between strains
- `panreactome_accumulation.tsv` -- reaction accumulation curve (pan-reactome growth with added genomes)
- `panreactome_dendrogram.png` -- hierarchical clustering of strains by metabolic similarity (with `--plot`)

The accumulation curve shows whether the pan-reactome is "open" (still growing) or "closed" (saturated), which has implications for how representative your strain collection is.

### Phenotype validation

Model predictions are only as trustworthy as their agreement with experimental data. The `validate` command compares in silico growth predictions against experimental growth phenotypes (e.g., Biolog plates, carbon source utilization assays) and computes standard classification metrics.

```bash
parmmosahn analyze validate \
  -M ./results/06-final_models/ \
  -p phenotypes.csv \
  -o ./analysis/validation.tsv
```

The phenotype file can be in **matrix format** (rows = carbon sources, columns = strains, values = 1/0) or **long format** (columns: `strain`, `carbon_source`, `growth`). Carbon sources should be specified as exchange reaction IDs (e.g., `EX_glc__D_e0`).

**Output metrics:** accuracy, sensitivity, specificity, and Matthews Correlation Coefficient (MCC). Per-strain per-substrate results are written to the TSV for detailed inspection.

### FBA summary

Run Flux Balance Analysis on all models to compare baseline growth rates and active exchange reactions:

```bash
parmmosahn analyze fba \
  -M ./results/06-final_models/ \
  -o ./analysis/fba_summary.tsv \
  -m medium.csv
```

### Auxotrophy screening

Systematically knock out each medium component to identify predicted auxotrophies -- strains that cannot grow without a specific nutrient. This reveals metabolic dependencies that may reflect ecological niche adaptation or gene loss events.

```bash
parmmosahn analyze auxotrophy \
  -M ./results/06-final_models/ \
  -o ./analysis/auxotrophy.tsv \
  -m medium.csv
```

### Gene essentiality prediction

Perform single-gene deletion FBA to predict which genes are essential for growth. This can be validated against experimental transposon library (Tn-seq) data and highlights potential drug targets in pathogens.

```bash
parmmosahn analyze essentiality \
  -M ./results/06-final_models/ \
  -o ./analysis/essentiality.tsv
```

### Reaction heatmap

Build a binary presence/absence matrix of all reactions across strains, optionally with a clustered heatmap visualization:

```bash
parmmosahn analyze heatmap \
  -M ./results/06-final_models/ \
  -o ./analysis/heatmap.tsv \
  --plot
```

## Modular usage

Each pipeline step can be run independently:

```bash
parmmosahn pan -g ./genomes/ -l my_label -t 8     # Steps 1-2 only
parmmosahn strains -g ./genomes/ -f soft.faa \
  -r rxns.tbl -m medium.csv                        # Step 3 only
parmmosahn evaluate -M ./models/ -o ./reports/     # Memote evaluation
parmmosahn annotate -M ./models/ -o ./annotated/   # FAIR annotation
parmmosahn gather -M ./models/ -o ./final/         # Format conversion
parmmosahn curate -s model.xml -c out/ -d curation.xlsx  # Apply curations
parmmosahn polish -M ./models/ -o ./polished/      # ModelPolisher
```

Run `parmmosahn --help` for a full list of commands and options.

## Configuration

PaRMMoSaHN supports YAML configuration files as an alternative to CLI flags
for the `run`, `pan`, and `project` commands. Generate a template:

```bash
parmmosahn init-config -o my_config.yml
```

The generated template documents every config-loadable key. Then use it:

```bash
# Full pipeline from a config file (all required args may come from config)
parmmosahn run --config my_config.yml

# Or mix CLI overrides with config defaults
parmmosahn run --config my_config.yml -o ./results/ -t 24
```

CLI arguments always override config file values. All options can also be set
via environment variables prefixed with `PARMMOSAHN_` (e.g.,
`PARMMOSAHN_THREADS=32`).

## ModelPolisher

ModelPolisher (v2.1-beta) enriches SBML models with BiGG database annotations
and standardised identifiers. It runs inside a container, so it requires
Docker, Podman, or Apptainer on the host.

**ModelPolisher is OFF by default**, both because of the container dependency
and because the bundled beta version of ModelPolisher uses a fragile network
fetch at startup. Enable it explicitly with `--polish`, optionally picking an
engine with `-e`:

```bash
# Default container engine (docker)
parmmosahn run ... --polish

# Podman (rootless, HPC-friendly)
parmmosahn run ... --polish -e podman

# Apptainer/Singularity (HPC clusters)
parmmosahn run ... --polish -e apptainer
```

If no container engine is detected at runtime, Step 5 is silently skipped and
the pipeline continues with the pre-polish models.

> **Note:** ModelPolisher v2.1 (stable release) has a known bug with a broken
> URL pattern regex in the DataONE namespace and crashes at startup. PaRMMoSaHN
> bundles the v2.1-beta which works correctly. SBML headers are temporarily
> downgraded from L3V2 to L3V1 for beta compatibility.

## Known Limitations

- **Reaction pre-filtering does not evaluate GPR rules.** If a reaction requires multiple gene subunits (AND rule), it may be included even if only one subunit has a homolog. The downstream `gapseq draft` step partially mitigates this.
- **Single medium for all strains.** Gap-filling uses one medium specification. Strains from different niches may need different media. Use `parmmosahn project` to re-derive models with alternative media.
- **Generic biomass composition.** Biomass reactions use gapseq's default gram-positive or gram-negative templates rather than species-specific composition. This is consistent with other automated tools (CarveMe, Bactabolize).
- **Soft-core threshold (2/N).** The default threshold excludes genes present in only one genome. For very small collections (N < 5), consider using `--add-unique` or adjusting `--soft-core`.

## Citation

If you use PaRMMoSaHN in your research, please cite:

```bibtex
@software{dewin2026parmmosahn,
  author = {De Win, Robbe and De Vrieze, Lucas},
  title = {PaRMMoSaHN: Pangenome Reference-based Metabolic Modelling by Saving Homolog Networks},
  version = {0.2.1},
  year = {2026},
  url = {https://github.com/robbedewin/PaRMMoSaHN}
}
```

## License

This project is licensed under the MIT License -- see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work was developed as part of a Master's thesis project at KU Leuven in the laboratory of Prof. Masschelein, in collaboration with the VIB-KU Leuven Center for Microbiology.
