# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [0.2.0] - 2026-04-14

### Added
- `parmmosahn analyze panreactome`: classify reactions as core/accessory/unique
  across a model collection; outputs Jaccard similarity matrix, accumulation
  curve (100 random orderings), per-reaction category TSVs, and optional
  Ward-linkage dendrogram PNG
- `parmmosahn analyze validate`: phenotype validation against experimental
  growth data; outputs per-strain confusion matrices, MCC, accuracy,
  sensitivity, and specificity
- `parmmosahn analyze essentiality`: single-gene deletion FBA for gene
  essentiality prediction; optionally compare against experimental Tn-seq /
  Keio knockout datasets
- `parmmosahn analyze fba`: FBA summary across a model collection (growth
  rates, active exchange reactions)
- `parmmosahn analyze auxotrophy`: screen all models for auxotrophies across
  amino acids and cofactors
- `parmmosahn analyze heatmap`: reaction presence/absence heatmap across
  strain models (requires `viz` extra)
- `--add-unique` flag on `run` and `pan`: optionally include unique/cloud gene
  families in the pan-model (default: off; soft-core only)
- Per-model biological statistics in `pipeline_summary.json`: reaction count,
  gene count, metabolite count per strain; core and accessory reaction totals
- Input genome validation: warns if any genome has fewer than 100 translated
  CDS features, or if fewer than 3 genomes are provided
- YAML config file support: persist pipeline parameters via `parmmosahn
  init-config` and `--config <file>`
- Step banners and indexed per-model progress indicators (`[i/n]`) in logging

### Fixed
- Bakta gene ID regex now handles both `G_gnl_Bakta_…` and `gnl_Bakta_…`
  formats; the annotation step was silently producing 0 NCBI protein
  accessions for all models built from PPanGGOLiN + Bakta-annotated genomes
- ModelPolisher binary detection now correctly matches the command actually
  invoked per engine (Docker vs Podman)
- SBML model gathering now uses recursive search so models nested in
  subdirectories are found correctly
- Model-identifier regex updated for compatibility with current gapseq output

### Changed
- `run` command writes `pipeline_summary.json` incrementally after each step;
  a crash mid-run now preserves timing and status for all completed steps
- Full INFO-level `run.log` appended to the output directory during `run` and
  `project`; all per-strain progress is captured even when the terminal session
  is lost
- `project` command now writes `pipeline_summary.json` on completion, making
  it the documented crash-recovery path for interrupted `run` executions
- `apply_curations` refactored into `_apply_duplicate_decisions` and
  `_apply_formula_fixes`; both the `curated_imbalances` sheet name (preferred)
  and the legacy `corrected_balances` sheet name are now accepted
- BLAST table filtering migrated to Polars for lower peak memory usage
- `run_command` uses threading for concurrent stdout/stderr capture, preventing
  deadlocks on large tool outputs
- `run_memote_snapshot` returns a boolean success flag for use in
  checkpointing logic

### Package
- `MANIFEST.in` and `pyproject.toml` `package-data` now list runtime data
  files explicitly; wildcard globs are removed to prevent large external
  reference files from being bundled in the wheel or sdist
- Thesis analysis scripts and genome data moved to the companion repository
  [PaRMMoSaHN-thesis](https://github.com/robbedewin/PaRMMoSaHN-thesis)


## [0.1.0] - 2026-03-11

### Added
- Initial release and core package structure
- Basic CLI interface for end-to-end execution
- Pangenome construction workflow (via PPanGGOLiN)
- Pan-model building capabilities (via gapseq)
- Strain model derivation and homolog network mapping
- "Human-in-the-loop" automatic curation tools (Excel-based)
- Docker/Podman integration for automated ModelPolisher execution
- Results gathering and formatting functionality
- Test suite (pytest), type checking (mypy), and code formatting (black/isort)
- CI/CD pipeline with GitHub Actions

[unreleased]: https://github.com/robbedewin/PaRMMoSaHN/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/robbedewin/PaRMMoSaHN/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/robbedewin/PaRMMoSaHN/releases/tag/v0.1.0
