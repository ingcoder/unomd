# EasyMD

EasyMD is a Python package that simplifies molecular dynamics simulations, making them accessible to both beginners and experts. It provides an automated, easy-to-use interface for running protein-ligand simulations using OpenMM as the backend.

![Demo Gif](assets/demo.gif)

## Who is it for?

- **Computational chemists** who want to streamline their MD workflow
- **Structural biologists** studying protein-ligand interactions
- **Drug discovery researchers** analyzing binding dynamics
- **Students and researchers** learning molecular dynamics
- **Anyone** who wants to run MD simulations without dealing with complex setup

## Key Features

- **Automated Setup**: From structure preparation to production runs
- **Integrated Force Fields**: AMBER14 for proteins, OpenFF 2.0.0 for small molecules
- **Flexible Configuration**: Easy to customize simulation parameters
- **Progress Monitoring**: Real-time updates on simulation progress
- **Analysis Tools**: Built-in tools for RMSD, RMSF calculations


## Tutorial

Follow our tutorial to learn how to use easy-md in quickrun mode. 
[Run the interactive tutorial in Google Colab](https://colab.research.google.com/drive/1H7IQ7mrGBOpuUN4-5XS8Vh09pwK3PuwL?usp=sharing)

## Prerequisites

First, install mamba, a fast package manager. You can use conda, but its extremely slow for some packages:

```bash
# If you have conda installed:
conda install mamba -n base -c conda-forge

# Or install mambaforge (standalone):
# For macOS/Linux:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

# For Windows:
# Download Mambaforge from:
# https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Windows-x86_64.exe
```

## Installation

### Using mamba (Recommended)

```bash
git clone https://github.com/ingcoder/easy-md.git
cd easy-md
mamba env create -f environment.yaml
mamba activate md_env
pip install -e 
```

## Quick Start

```python
from easy_md.main.quickrun import quickrun

# Run a simple protein-ligand simulation
quickrun(
    protein_file="path/to/protein.pdb",
    ligand_file="path/to/ligand.sdf",
    nsteps=1000
)
```

## Requirements

Core dependencies are automatically managed through pyproject.toml:

- Python >=3.9
- OpenMM >=7.7.0
- OpenFF-Toolkit >=0.11.0
- MDAnalysis >=2.4.0

Additional development dependencies are available in the conda environment.yaml file.

## Project Structure

```
my_project/
├── config/
│   └── simulation_config.yaml  # Simulation parameters
├── structures/
│   ├── protein.pdb            # Input protein structure
│   └── ligand.sdf            # Input ligand structure
└── output/                    # Simulation outputs
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

We welcome contributions! Please feel free to submit a Pull Request.

## Citation

If you use EasyMD in your research, please cite it as:

```bibtex
@software{easymd2025,
  author       = {Ingrid Barbosa-Farias, SimAtomic},
  title        = {EasyMD: A Python Package for Simplified Molecular Dynamics Simulations},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/ingcoder/easy-md}
}
```
