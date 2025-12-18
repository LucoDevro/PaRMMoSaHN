# PaRMMoSaHN

**Pa**ngenome-based **R**econstruction of **M**etabolic **Mo**dels of **S**trains **a**nd **H**omologous **N**etworks

[![CI](https://github.com/robbedewin/PaRMMoSaHN/workflows/CI/badge.svg)](https://github.com/robbedewin/PaRMMoSaHN/actions)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> 🚧 **Status: Under Active Development** - Documentation and features are being actively implemented. See [docs/TODO.md](docs/TODO.md) for tracking progress.

## Overview

PaRMMoSaHN is a Python library for pangenome-based reconstruction of strain-specific metabolic models from bacterial genomes. It enables researchers to construct accurate metabolic models by leveraging pangenome information to capture strain-level diversity and generate homologous network representations.

## Features

<!-- TODO: Update with actual implemented capabilities as they're completed -->

- [ ] Construct pangenome-based metabolic models (in progress)
- [ ] Derive strain-specific models from pangenome references (planned)
- [ ] Automated and manual curation workflows (planned)
- [ ] Integration with standard metabolic modeling frameworks (planned)

## Installation

### From Source

Clone the repository and install in development mode:

```bash
git clone https://github.com/robbedewin/PaRMMoSaHN.git
cd PaRMMoSaHN
pip install -e .
```

### Dependencies

PaRMMoSaHN requires Python 3.9 or higher. Key dependencies include:
- numpy
- pandas
- scikit-learn
- matplotlib

For a complete list, see [`requirements.txt`](requirements.txt).

## Quick Start

### Workflow

1. **Construct Pangenome**: Build a pangenome from genome sequences
2. **Build Pan-Model**: Create a pangenome-wide metabolic model
3. **Derive Strain Models**: Generate strain-specific models
4. **Curate Models**: Apply manual and automatic curation
5. **Analyze Results**: Gather and visualize results

### Example Usage

<!-- TODO: Replace with real usage examples once core functionality is implemented -->

```python
# Example placeholder - will be updated with actual API
from package.core import PangenomeBuilder, ModelDeriver

# TODO: Add real workflow example here
```

For command-line interface:

```bash
parmmosahn --help  # TODO: Document CLI commands once implemented
```

See the [notebooks/](notebooks/) directory for detailed workflow examples (in progress).

## Contributing

Contributions are welcome! Please follow these steps to contribute:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature/YourFeature`).
3. Make your changes and commit them (`git commit -m 'Add some feature'`).
4. Push to the branch (`git push origin feature/YourFeature`).
5. Open a pull request.

## Citation

If you use PaRMMoSaHN in your research, please cite:

```bibtex
@software{dewin2025parmmosahn,
  author = {De Win, Robbe},
  title = {PaRMMoSaHN: Pangenome Reference-based Metabolic Modelling by Saving Homolog Networks},
  version = {0.1.0},
  year = {2025},
  url = {https://github.com/robbedewin/PaRMMoSaHN}
}
```

See [`CITATION.cff`](CITATION.cff) for more citation formats.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments
This work was developed as part of a thesis project at KU Leuven in the laboratory of Prof. Masschelein, in collaboration with the VIB-KU Leuven Center for Microbiology.

## Contact

For questions or issues, please [open an issue](https://github.com/robbedewin/PaRMMoSaHN/issues) on GitHub.
