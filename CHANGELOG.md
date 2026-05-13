# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [0.2.1] - 2026-05-13

### Added
- `pipeline_summary.json` `provenance` block: SHA-256 of the medium CSV
  and curation database, plus the soft-core threshold and full PPanGGOLiN
  configuration actually used. Lets a reviewer reconstruct any v0.2.1 run
  from the summary alone.
- `pipeline_summary.json` `run_environment` block: host CPU count,
  total and available RAM, platform string, Python version.
- `01-pangenome/pangenome_meta.json` sidecar written during pangenome
  construction; records the parameters that produced the soft-core
  proteome for provenance recovery by `write_pipeline_summary`.
- `04-curated_models/curation_application_report.tsv` listing per-row
  hit counts across strains; rows that fired on 0 strains emit a
  WARNING. Addresses the silent-no-op behaviour of stale curation-DB
  sheets documented during thesis validation.
- `parmmosahn doctor` "Host resources" block reporting CPU/RAM/OS, a
  recommended memote-worker count, and a low-RAM warning below 8 GB.
- `parmmosahn evaluate --workers N / --skip-existing` flags for parallel
  memote evaluation with per-strain checkpointing; runtime warning when
  `--workers N` exceeds the ~4 GB/worker memote memory budget — surfaces
  the WSL2 `BrokenProcessPool` footgun that previously required reading
  the docs to discover.
- `apply_curations` now accepts an optional `stats` accumulator dict so
  callers can aggregate per-row hit counts across strains without
  re-walking the SBML files. Public signature is otherwise unchanged.
- Helpers `parmmosahn.utils.sha256_file` and
  `parmmosahn.utils.get_system_resources`.

### Changed
- ModelPolisher (Step 5) is now **opt-in** rather than opt-out. The CLI
  flag changed from `--no-polish` (polishing ON by default) to `--polish`
  (polishing OFF by default). Rationale: the bundled ModelPolisher 2.1-beta
  performs a fragile MIRIAM/identifiers.org network fetch at startup that
  intermittently fails, the container engine is an extra dependency, and
  running polishing on PaRMMoSaHN but not on comparator pipelines (CarveMe,
  gapseq) inflates memote's annotation sub-score for the wrong reason. This
  is the only user-visible flag change in v0.2.1 — anyone scripting
  `--no-polish` will need to remove that flag; anyone who wants polishing
  must now pass `--polish` explicitly.
- `derive_strain_models` now cancels queued strain jobs on `KeyboardInterrupt`
  rather than waiting for the ThreadPoolExecutor to drain naturally. On
  multi-hour pipeline runs this means `Ctrl-C` returns a prompt within
  seconds instead of after the next strain finishes.
- `EXPECTED_STRAINS` (tests) updated to four members of the Neal et al.
  2024 *B. subtilis* cohort: `GCF_000009045.1`, `GCF_000155325.1`,
  `GCF_000209795.2`, `GCF_000227485.1`. The previous list mixed a stray
  `GCA_000227485.1` (GenBank prefix; every other entry was RefSeq) with
  a `GCF_000146565.1` that was not present in any local artefact.
- `cli.run` and `cli.project` parameter dictionaries now include
  `curation_db`, `add_unique`, `diamond_bits`, `gapseq_bits`, `biomass`,
  and `parallel_strains` so the summary captures the complete invocation.

### Tests
- Added 11 unit tests covering the new provenance helpers, summary
  fields, curation-applicability report, and the doctor host-resources
  block.
- Added integration tests under the `biological` mark that verify
  `pipeline_summary.json` provenance and `pangenome_meta.json` sidecar
  contents on a real 4-genome run.

### Notes
- v0.2.1 contains no changes to reaction stoichiometry, gap-fill logic,
  curation rules, or annotation enrichment. Every SBML produced by
  v0.2.0 is reproduced byte-for-byte by v0.2.1, so all validation
  results in the companion thesis repository remain valid.


## [0.2.0] - 2026-04-15

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
- `parmmosahn analyze validate` no longer leaks carbon sources between
  phenotype tests: a per-strain baseline medium is snapshotted before the
  first test and restored for each subsequent one. Previously, accuracy/MCC
  numbers were silently inflated because carbon sources from earlier tests
  remained enabled on `model.medium`
- `parmmosahn analyze validate` growth-column parser now accepts bool,
  int/float, and common string forms (`True`/`False`, `yes`/`no`, `1`/`0`);
  previously, string booleans raised `ValueError` despite being documented
  as valid input

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

### Removed
- Stale `errors.md` troubleshooting log (contained local hostname/path
  references and was superseded by the in-output `run.log` and `errors.log`
  files written by the `run` and `project` commands)

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

[unreleased]: https://github.com/robbedewin/PaRMMoSaHN/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/robbedewin/PaRMMoSaHN/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/robbedewin/PaRMMoSaHN/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/robbedewin/PaRMMoSaHN/releases/tag/v0.1.0
