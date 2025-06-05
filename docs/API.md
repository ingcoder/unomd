# EasyMD API Documentation

This document provides detailed information about the key functions and classes in EasyMD.

## Main Module Functions

### quickrun

```python
def quickrun(protein_file: str, ligand_file: str = None, nsteps: int = 1000) -> None
```

Runs a complete molecular dynamics simulation workflow with minimal configuration.

**Parameters:**
- `protein_file` (str): Path to the input protein PDB file
- `ligand_file` (str, optional): Path to the ligand SDF file
- `nsteps` (int, optional): Number of simulation steps (default: 1000)

**Example:**
```python
from easy_md.main.quickrun import quickrun

quickrun(
    protein_file="protein.pdb",
    ligand_file="ligand.sdf",
    nsteps=1000
)
```

### run_solvation.add_water

```python
def add_water(config: dict) -> None
```

Solvates the molecular system by adding water molecules and ions.

**Parameters:**
- `config` (dict): Configuration dictionary containing simulation parameters

**Key Configuration Parameters:**
- `solv_box_buffer`: Buffer size around the solute (in Å)
- `solv_ionic_strength`: Ionic concentration (in M)
- `solv_positive_ion`: Type of positive ion (Na+, K+, etc.)
- `solv_negative_ion`: Type of negative ion (Cl-, etc.)

### run_forcefield_parameterization.main

```python
def main(config: dict, print_detailed_info: bool = False) -> None
```

Applies force field parameters to the molecular system.

**Parameters:**
- `config` (dict): Configuration dictionary
- `print_detailed_info` (bool): Whether to print detailed force field information

### run_energy_minimization.main

```python
def main(config: dict) -> None
```

Performs energy minimization of the molecular system.

**Parameters:**
- `config` (dict): Configuration dictionary containing simulation parameters

### run_simulation.main

```python
def main(config: dict = None, starting_state_path: str = None, starting_checkpoint_path: str = None, equilibration_only: bool = False) -> None
```

Runs the production molecular dynamics simulation.

**Parameters:**
- `config` (dict): Configuration dictionary
- `starting_state_path` (str, optional): Path to a starting state file
- `starting_checkpoint_path` (str, optional): Path to a checkpoint file
- `equilibration_only` (bool): Whether to run only the equilibration phase

## Utility Functions

### config.create_config

```python
def create_config(
    protein_file: str = None,
    ligand_file: str = None,
    project_dir: str = None,
    output_dir: str = None,
    config_dir: str = None,
    save_config_as: str = "simulation_config.yaml",
    **params
) -> Dict[str, Any]
```

Creates a configuration dictionary with default and user-provided settings.

**Parameters:**
- `protein_file` (str): Path to the protein PDB file
- `ligand_file` (str, optional): Path to the ligand file
- `project_dir` (str, optional): Project directory path
- `output_dir` (str, optional): Output directory path
- `config_dir` (str, optional): Configuration directory path
- `save_config_as` (str): Name of the configuration file
- `**params`: Additional configuration parameters

**Returns:**
- Dict[str, Any]: Configuration dictionary

### simulation_util.setup_barostat

```python
def setup_barostat(temperature: float, pressure: float, barostat_frequency: int, use_anisotropic: bool = False) -> Union[MonteCarloBarostat, MonteCarloAnisotropicBarostat]
```

Configures and returns a barostat for pressure control.

**Parameters:**
- `temperature` (float): System temperature in Kelvin
- `pressure` (float): Target pressure in atmospheres
- `barostat_frequency` (int): Frequency of barostat updates
- `use_anisotropic` (bool): Whether to use anisotropic pressure coupling

**Returns:**
- OpenMM Barostat object

## Configuration Parameters

### Paths
- `path_base`: Base project directory
- `path_protein`: Input protein structure file
- `path_ligand`: Input ligand structure file
- `path_protein_solvated`: Solvated system output
- `path_openff_topology`: OpenFF topology file
- `path_openmm_system`: OpenMM system file
- `path_emin_structure`: Energy minimized structure
- `path_md_trajectory`: Trajectory output file
- `path_md_checkpoint`: Checkpoint file

### Force Fields
- `ff_small_molecule_openff`: Small molecule force field
- `ff_protein_openff`: Protein force field
- `ff_protein`: Alternative protein force field
- `ff_water`: Water model

### Integrator Settings
- `integrator_temperature`: System temperature (K)
- `integrator_friction`: Friction coefficient (1/ps)
- `integrator_timestep`: Integration time step (ps)

### Equilibration Parameters
- `total_steps`: Total number of steps
- `save_interval`: Frequency of saving coordinates
- `pressure_atm`: Target pressure
- `barostat_freq`: Barostat update frequency

### Platform Settings
- `platform_name`: Computation platform (CUDA, OpenCL, CPU)
- `platform_properties`: Platform-specific settings 