# Contributing to PaRMMoSaHN

Thank you for your interest in contributing to PaRMMoSaHN! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/PaRMMoSaHN.git
   cd PaRMMoSaHN
   ```
3. Create a virtual environment and install dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -e ".[dev]"
   ```
4. Install pre-commit hooks (this will run linters and formatters before each commit):
   ```bash
   pre-commit install
   ```

## Development Workflow

### Branching Strategy

- `main` - Stable, production-ready code
- `develop` - Integration branch for features (if needed)
- `feature/*` - New features
- `bugfix/*` - Bug fixes
- `hotfix/*` - Urgent fixes for production

### Making Changes

1. Create a new branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes following the code style guidelines

3. Run tests and linting:
   ```bash
   pytest tests/
   black src/ tests/
   isort src/ tests/
   mypy src/
   ```

4. Commit your changes with descriptive messages:
   ```bash
   git commit -m "Add feature: brief description"
   ```

5. Push to your fork and submit a pull request

## Code Style

- Follow PEP 8 guidelines
- Use type hints for function signatures
- Write docstrings for all public functions and classes
- Maximum line length: 88 characters (Black default)
- Use `black` for code formatting
- Use `isort` for import sorting

## Testing

- Write tests for all new features and bug fixes
- Ensure all tests pass before submitting a PR
- Aim for high test coverage
- Use pytest for writing tests

### Biological integration tests

The `pytest -m biological` suite validates real pipeline outputs against
biological invariants (FBA feasibility, reaction-count ranges, compartment
layout, regression checks for past bugs). It is **opt-in** and excluded
from the default `pytest` invocation because it requires:

1. A `test_run/` directory at the repo root with pre-computed pangenome +
   pan-model + strain-model artefacts.
2. External tools available on `$PATH`: PPanGGOLiN, gapseq, DIAMOND.

The strain set is pinned in [`tests/conftest.py`](tests/conftest.py): four
RefSeq members of the Neal et al. 2024 *B. subtilis* cohort
(`GCF_000009045.1`, `GCF_000155325.1`, `GCF_000209795.2`, `GCF_000227485.1`).
Each test gracefully skips for strains that are missing from `test_run/`,
so partial coverage is acceptable.

#### Bootstrapping `test_run/` from a fresh run

```bash
# 1. Stage the four genomes + medium CSV
mkdir -p test_run/genomes
for s in GCF_000009045.1 GCF_000155325.1 GCF_000209795.2 GCF_000227485.1; do
  ln -sf /path/to/$s.gbff test_run/genomes/$s.gbff
done
ln -sf /path/to/M9_glc_aerobic.csv test_run/M9_glc_aerobic.csv

# 2. Run the pipeline through Step 3 (~50-90 min on 4 small genomes)
parmmosahn run \
  -g test_run/genomes -o test_run -l test_run \
  -m test_run/M9_glc_aerobic.csv \
  -t 8 --parallel-strains 2 --pause-for-curation

# 3. Run the suite
pytest -m biological
```

The `--pause-for-curation` flag halts after Step 3 (Steps 4–7 are not
exercised by the biological tests and add ~30 min of memote runtime).
`test_run/` is in `.gitignore` and should not be committed.

For long pipeline runs, wrapping the `parmmosahn run` invocation in
`tmux new -s parmmosahn_testrun '...'` is recommended so the run
survives terminal detach.

## Documentation

- Update the README.md if you change functionality
- Update docstrings for any modified functions/classes
- Add entries to CHANGELOG.md for notable changes
- Update type hints if function signatures change

## Pull Request Process

1. Update the CHANGELOG.md with details of changes
2. Ensure all tests pass and pre-commit hooks succeed
3. Update documentation as needed
4. Fill out the pull request template completely
5. Request review from maintainers
6. Address any review comments

## Code Review

- Be respectful and constructive in reviews
- Focus on the code, not the person
- Explain your reasoning for requested changes
- Be open to feedback on your own code

## Reporting Bugs

Use the bug report template when creating issues. Include:
- Clear description of the bug
- Steps to reproduce
- Expected vs actual behavior
- Environment details (OS, Python version, etc.)
- Error messages/tracebacks

## Feature Requests

Use the feature request template. Describe:
- The problem you're trying to solve
- Your proposed solution
- Alternative solutions considered
- Whether you're willing to implement it

## Questions?

Feel free to open an issue with your question or reach out to the maintainers.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
