# PaRMMoSaHN

**Pa**ngenome **R**eference-based **M**etabolic **Mo**delling by **Sa**ving **H**omolog **N**etworks

[![CI](https://github.com/robbedewin/PaRMMoSaHN/workflows/CI/badge.svg)](https://github.com/robbedewin/PaRMMoSaHN/actions)
[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**PaRMMoSaHN** is a Python-based pipeline that bridges partitioned pangenome graphs with high-fidelity metabolic reconstruction.

Instead of building metabolic models for hundreds of strains from scratch (which is computationally expensive and prone to inconsistent gap-filling), PaRMMoSaHN builds a singular high-quality "pan-model" and rapidly derives strain-specific models by mapping strains to the pan-model via homology networks.

## Features

- **Pangenome Construction:** Automated soft-core genome extraction using [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN).
- **Pan-Model Generation:** High-fidelity pangenome-scale metabolic network reconstruction via [gapseq](https://gapseq.readthedocs.io/).
- **Rapid Strain Derivation:** Maps individual strains to the pan-model using [DIAMOND](https://github.com/bbuchfink/diamond) and gap-fills missing reactions.
- **Automated Curation:** Handles duplicates, stoichiometry correction, and MEMOTE-driven curation via an intuitive Excel decision database using `cobrapy`.
- **FAIR Compliance:** Optionally polishes and annotates final SBML models using the Palsson lab's [ModelPolisher](https://github.com/draeger-lab/ModelPolisher) (supports Docker, Podman, and Apptainer for HPCs).
- **Multi-format Export:** Exports final validated models to SBML (.xml), JSON, and MATLAB formats.

## Installation

Because PaRMMoSaHN orchestrates several complex bioinformatics binaries, **you must use Conda/Mamba** to install the environment. Standard `pip` cannot install tools like `diamond` or `gapseq`.

### 1. Create the Conda Environment
Download the `environment.yml` file from this repository. This file strictly pins the compatible versions of all external dependencies. *(Note: We highly recommend using `mamba` or `micromamba` for faster installation).*

```bash
conda env create -f environment.yml
conda activate parmmosahn_env

```

### 2. Install PaRMMoSaHN

*(Currently in Private Beta)* Install the package directly from GitHub:

```bash
pip install "git+[https://github.com/robbedewin/PaRMMoSaHN.git@v0.1.0](https://github.com/robbedewin/PaRMMoSaHN.git@v0.1.0)"

```

## Quick Start

PaRMMoSaHN is primarily driven through its Command Line Interface (CLI). There are two main ways to run the pipeline, depending on whether you want to manually curate the draft models.

### Workflow A: Fully Automated Run (No Manual Curation)
To execute the pipeline from start to finish, automatically carrying forward draft models without manual intervention, simply run:

```bash
parmmosahn run \
  --genomes-dir /path/to/gbff_files/ \
  --output-dir ./parmmosahn_results/ \
  --label My_Strains \
  --medium medium_config.csv \
  --threads 14
```

### Required Inputs

1. **`--genomes-dir`**: A directory containing your strain genomes in GenBank format (`.gbff`).
2. **`--medium`**: A `.csv` file defining the growth medium for gap-filling (gapseq format).
3. **`--label`**: A string label for the pangenome/pan-model output files (e.g., 'clostridiales').


### Workflow B: Human-in-the-Loop Curation

For highly accurate models, PaRMMoSaHN allows you to pause the pipeline, manually evaluate the draft models using an Excel spreadsheet, and then resume the pipeline to polish and export them.

**Step 1: Generate Draft Models and Pause**

Use the `--stop-after-drafts` flag. The pipeline will build the pangenome, pan-model, and draft strains, and then stop.

```bash
parmmosahn run \
  --genomes-dir ./data \
  --output-dir ./results \
  --label My_Strains \
  --medium medium.csv \
  --stop-after-drafts

```

**Step 2: Generate the Evaluation Spreadsheet**

Use the `evaluate` commands to run MEMOTE tests on the draft models and extract the problems into a human-readable Excel template.

```bash
parmmosahn evaluate snapshot --model-dir ./results/03-strains --out-dir ./reports
parmmosahn evaluate extract --report-dir ./reports --model-dir ./results/03-strains

```

*Open the resulting `curation_template.xlsx` in Excel, make your curation decisions in the `curated_imbalances` and `duplicate_reactions` sheets, and save the file.*

**Step 3: Resume the Pipeline**

Provide your filled Excel database to the `resume` command to apply your curations, run ModelPolisher, and generate the final FAIR-compliant SBML models.

```bash
parmmosahn resume \
  --drafts-dir ./results/03-strains \
  --output-dir ./results \
  --curation-db curation_template.xlsx

```

## Advanced Configuration & HPC Usage

PaRMMoSaHN exposes advanced parameters for highly diverse pangenomes and High-Performance Computing (HPC) environments.

```bash
parmmosahn run ... \
  --cluster-id 0.80 \
  --cluster-cov 0.80 \
  --diamond-bits 150 \
  --gapseq-bits 150

```

### ModelPolisher on HPC Clusters

By default, Step 6 (ModelPolisher) attempts to use `docker compose`. If you are running on an HPC cluster where Docker is forbidden, you can switch the container engine to **Apptainer** or **Podman**, or skip the step entirely:

```bash
# Use Apptainer (Singularity) instead of Docker
parmmosahn run ... --container-engine apptainer

# Or skip the polishing step entirely
parmmosahn run ... --skip-polisher

```

## Modular Usage

You can also run specific parts of the pipeline individually if you only need a pangenome or a pan-model:

```bash
parmmosahn pangenome --genomes-dir ./data --label test --threads 8
parmmosahn panmodel --fasta softcore.faa --threads 8

```

Run `parmmosahn --help` for a full list of commands.

## Citation

If you use PaRMMoSaHN in your research, please cite:

```bibtex
@software{dewin2026parmmosahn,
  author = {De Win, Robbe and De Vrieze, Lucas},
  title = {PaRMMoSaHN: Pangenome Reference-based Metabolic Modelling by Saving Homolog Networks},
  version = {0.1.0},
  year = {2026},
  url = {https://github.com/robbedewin/PaRMMoSaHN}
}

```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work was developed as part of a Master's thesis project at KU Leuven in the laboratory of Prof. Masschelein, in collaboration with the VIB-KU Leuven Center for Microbiology.
