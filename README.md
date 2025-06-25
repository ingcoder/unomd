<!-- or percentage of the text‐column width -->
<p align="center">
  <img src="assets/feynman_quote_3.png"
       alt="Feynman quote"
       width=70%>
</p>

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
mamba activate easymd
pip install -e .
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

## Step-by-Step Approach
```python
from easy_md.utils.config import create_config
from easy_md.main import run_solvation, run_forcefield_parameterization, run_energy_minimization, run_simulation

config = create_config(
    protein_file="path/to/protein.pdb",
    ligand_file="path/to/ligand.sdf",
    
    # MD simulation settings. See "Simulation Paramters" section below for all options
    md_steps=1000,
    md_save_interval=10,
    
    # Platform settings
    platform_name="CPU",          # or CUDA
    platform_precision="mixed",   # or "single" or "double"
)
run_solvation.add_water(config=config)
run_forcefield_parameterization.main(config)
run_energy_minimization.main(config)
run_simulation.main(config, starting_state_path="path/to/state.xml")
# By default `run_simulation.main(config)` loads the energy-minimized state
# saved in `emin.xml`.  To resume a previous run instead, supply the path
# to its state file .xml starting_state_path="path/to/state.xml" or checkpoint file: starting_state_path="path/to/state.xml":
```


## Simulation Parameters
```yaml
# Energy Minimization Parameters
emin_heating_interval: 1           # Interval for heating during minimization
emin_heating_step: 300             # Heating step size
emin_steps: 10                     # Number of minimization steps
emin_target_temp: 300              # Target temperature for minimization (K)
emin_tolerance: 5                  # Energy tolerance (kJ/mol/nm)

# Force Field Parameters
ff_protein: "amber14-all.xml"                               # Protein force field
ff_protein_openff: "ff14sb_off_impropers_0.0.3.offxml"      # OpenFF protein parameters
ff_small_molecule_openff: "openff-2.0.0.offxml"             # Small molecule force field
ff_water: "amber14/tip3pfb.xml"                             # Water model

# Integrator Settings
integrator_friction: 1.0          # Langevin integrator friction coefficient (1/ps)
integrator_temperature: 300.0     # Simulation temperature (K)
integrator_timestep: 0.002        # Integration timestep (ps)

# Molecular Dynamics Settings
md_anisotropic: false             # Use anisotropic barostat
md_barostat_freq: 25              # Barostat update frequency
md_harmonic_restraint: true       # Apply harmonic restraints
md_load_state: true               # Load from previous state if available
md_npt: false                     # Run in NPT ensemble
md_pressure: 1.0                  # System pressure (atm)
md_restrained_residues: []        # List of residues to restrain
md_save_interval: 10              # Trajectory save interval
md_steps: 1000                    # Number of MD steps

# Monitoring Parameters
monitor_energy_threshold: 100.0    # Energy monitoring threshold
monitor_temp_threshold: 2.0        # Temperature monitoring threshold
monitor_window: 10                 # Monitoring window size

# Solvation Parameters
solv_box_buffer: 2.5              # Solvent box padding (Å)
solv_ionic_strength: 0.15         # Ionic strength (M)
solv_model: "tip3p"               # Water model
solv_negative_ion: "Cl-"          # Negative ion type
solv_pH: 7.0                      # System pH
solv_positive_ion: "Na+"          # Positive ion type

# Analysis Settings
rmsd_ligand_selection: "resname UNK"   # RMSD selection for ligand
rmsd_selection: "protein and name CA"  # RMSD selection for protein
rmsf_selection: "protein and name CA"  # RMSF selection

# Platform Configuration
platform_name: "GPU"               # Computation platform (GPU/CPU)
platform_precision: "mixed"        # Precision model

# Paths
path_protein: path/to/protein.pdb                # Required by user
path_ligand: path/to/ligand.sdf                  # Required if you want to simulate protein-ligand interaction

# Paths set automatically, unless provided by user

path_base: path/to/project_dir # # Defaults to the parent directory containing the protein structure files

path_amber_topology: path_base/output/amber_top.prmtop
path_emin_state: path_base/output/emin.xml
path_emin_structure: path_base/output/emin.pdb
path_md_checkpoint: path_base/output/md_checkpoint_id.chk
path_md_image: path_base/output/md_image_id.dcd                #! Final processed trajectory file with molecules re-centered in water box.
path_md_log: path_base/output/md_id.log
path_md_state: path_base/output/md_state_id.xml
path_md_trajectory: path_base/output/md_trajetory_id.dcd
path_openff_interchange: path_base/output/openff_interchange.pdb
path_openff_topology: path_base/output/openff_topology.json
path_openmm_system: path_base/output/openmm_system.xml
path_openmm_topology: path_base/output/openmm_topology.pkl
path_protein_solvated: path_base/output/protein_solvated.pdb
path_rmsd_ligand_output: path_base/output/rmsd_ligand.pkl
path_rmsd_output: path_base/output/rmsd.pkl
path_rmsf_output: path_base/output/rmsf.log
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
