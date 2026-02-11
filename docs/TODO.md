# Documentation TODO

This file tracks documentation tasks to complete before publication.

## Priority: Before v0.1.0 Release

### README.md
- [ ] Update "Features" with actual implemented capabilities
- [ ] Add real usage examples (replace placeholder code)
- [ ] Add workflow diagram or architecture overview
- [ ] Update installation instructions once dependencies are finalized
- [ ] Add troubleshooting section if needed
- [ ] Add performance/system requirements
- [ ] Update example notebooks links once workflows are complete

### CITATION.cff
- [ ] Add your real email address (replace `your.email@example.com`)
- [ ] Add your ORCID (register at https://orcid.org if needed)
- [ ] Update date-released when v0.1.0 is ready
- [ ] Add DOI once Zenodo release is created
- [ ] Add Lucas De Vrieze as co-author

### CONTRIBUTING.md
- [ ] Update development workflow once stable
- [ ] Add specific code style examples from actual codebase
- [ ] Document how to run full test suite once implemented
- [ ] Add architecture/design docs reference if created

### Code Documentation
- [ ] Add docstrings to all public functions in src/package/
- [ ] Create API documentation (neccessary?)
- [ ] Document CLI commands once implemented
- [ ] Add type hints throughout codebase

### Notebooks
- [ ] Complete notebooks/01-construct_pangenome.ipynb
- [ ] Complete notebooks/02-build_pan_model.ipynb
- [ ] Complete notebooks/03-derive_strain_models.ipynb
- [ ] Complete notebooks/04-manual_curations.ipynb
- [ ] Complete notebooks/05-automatic_curations.ipynb
- [ ] Complete notebooks/06-gather_results.ipynb
- [ ] Add narrative explanations in each notebook
- [ ] Ensure notebooks run end-to-end without errors

## Nice to Have (Post v0.1.0)

### Documentation Site
- [ ] Set up Sphinx or MkDocs for full documentation
- [ ] Add tutorials and how-to guides
- [ ] Add API reference
- [ ] Add theoretical background/methods section

### Media
- [ ] Create logo or visual identity
- [ ] Add workflow diagrams
- [ ] Create demo video or GIF
- [ ] Add result visualizations

### Community
- [ ] Write blog post or preprint
- [ ] Create discussion forum or Slack/Discord
- [ ] Add FAQ section
- [ ] Create roadmap document

### Zenodo/Publication
- [ ] Create Zenodo account and link repository
- [ ] Create v0.1.0 release on GitHub
- [ ] Get DOI from Zenodo
- [ ] Update all references to include DOI
- [ ] Consider submitting to JOSS or similar

## Content Templates

### For README Usage Section
```python
# Once implemented, replace with:
from package.core import PangenomeBuilder
from package.modeling import ModelDeriver

# Build pangenome
builder = PangenomeBuilder(genome_dir="data/00_fasta/pan/")
pangenome = builder.construct()

# Create pan-model
model = builder.create_pan_model(pangenome)

# Derive strain-specific model
strain_model = ModelDeriver.derive(model, strain_id="GCF_000016505")
```

### For CLI Documentation
```bash
# Once implemented:
parmmosahn build-pangenome --input data/00_fasta/pan/ --output data/01_pangenome/
parmmosahn create-model --pangenome data/01_pangenome/ --output models/
parmmosahn derive-strain --model models/pan_model.xml --strain GCF_000016505
```

## Notes
- Keep CHANGELOG.md updated as features are implemented
- Update version numbers in pyproject.toml for each release
- Test installation instructions on fresh environment before release
- Run all pre-commit hooks before committing documentation updates
