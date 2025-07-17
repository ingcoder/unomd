# Getting Started with UnoMD

This tutorial will guide you through common use cases of UnoMD, from basic protein-ligand simulations to more advanced scenarios.

## Basic Protein-Ligand Simulation

### 1. Quick Start

The simplest way to run a simulation is using the `quickrun` function:

```python
from unomd.main.quickrun import quickrun

quickrun(
    protein_file="protein.pdb",
    ligand_file="ligand.sdf",
    nsteps=1000
)
```

This will automatically:

1. Solvate your system
2. Apply force field parameters
3. Perform energy minimization
4. Run a short simulation

### 2. Step-by-Step Approach

For more control over the process, you can run each step individually:

```python
from unomd.main import run_solvation, run_forcefield_parameterization
from unomd.main import run_energy_minimization, run_simulation
from unomd.utils.config import create_config

# Create configuration
config = create_config(
    protein_file="protein.pdb",
    ligand_file="ligand.sdf",
    project_dir="my_simulation"
)

# Step 1: Solvate the system
run_solvation.add_water(config=config)

# Step 2: Apply force fields
run_forcefield_parameterization.main(config)

# Step 3: Minimize energy
run_energy_minimization.main(config)

# Step 4: Run simulation
run_simulation.main(config)
```

## Advanced Usage

### 1. Customizing Simulation Parameters

Create a custom configuration file (`config.yaml`):

```yaml
paths:
  base_folder: "my_simulation"
  ligand: "ligand.sdf"
  solvated_protein: "protein_solvated.pdb"

integrator:
  temperature_kelvin: 310 # Body temperature
  friction_coeff_ps: 1
  time_step_ps: 0.002

equilibration:
  total_steps: 5000000 # 10 ns
  save_interval: 50000 # Save every 100 ps
  pressure_atm: 1.0
  barostat_freq: 25
```

Load and use the custom configuration:

```python
import yaml
from unomd.utils.config import create_config

# Load custom settings
with open("config.yaml") as f:
    custom_settings = yaml.safe_load(f)

# Create config with custom settings
config = create_config(**custom_settings)
```

### 2. Continuing from a Checkpoint

```python
from unomd.main import run_simulation

# Continue from a checkpoint
run_simulation.main(
    config=config,
    starting_checkpoint_path="output/md_checkpoint_0.chk"
)
```

### 3. Running NPT Equilibration

For membrane proteins or when volume equilibration is important:

```python
config = create_config(
    protein_file="protein.pdb",
    ligand_file="ligand.sdf",
    md_npt=True,
    md_pressure=1.0,  # atm
    md_barostat_freq=25
)

run_simulation.main(config)
```

### 4. Using Position Restraints

To restrain specific residues during simulation:

```python
config = create_config(
    protein_file="protein.pdb",
    ligand_file="ligand.sdf",
    md_harmonic_restraint=True,
    md_restrained_residues=[1, 2, 3, 4, 5]  # Residue indices to restrain
)

run_simulation.main(config)
```

## Tips and Best Practices

1. **System Preparation**

   - Always check your input structures
   - Ensure proper protonation states
   - Remove any unwanted molecules/ions

2. **Simulation Length**

   - Start with short test runs (1000 steps)
   - Increase length for production runs
   - Monitor energy and temperature convergence

3. **File Management**

   - Use descriptive file names
   - Keep input files organized
   - Save checkpoints regularly

4. **Performance Optimization**
   - Use CUDA platform when available
   - Adjust save intervals based on needs
   - Consider periodic boundary conditions

## Troubleshooting

### Common Issues

1. **System Crashes**

   ```python
   # Reduce time step
   config = create_config(
       protein_file="protein.pdb",
       integrator_timestep=0.001  # Reduced from 0.002
   )
   ```

2. **Energy Instability**

   ```python
   # Add more minimization steps
   run_energy_minimization.main(config, max_iterations=5000)
   ```

3. **Memory Issues**
   ```python
   # Reduce trajectory saving frequency
   config = create_config(
       protein_file="protein.pdb",
       md_save_interval=5000  # Save less frequently
   )
   ```

### Getting Help

- Check the [API Documentation](API.md)
- Look for similar issues in the repository
- Contact the developers with:
  - Your configuration file
  - Error messages
  - System details
