"""This script performs simulation of a molecular system using OpenMM."""

# Standard library imports
from collections import deque
from pathlib import Path

# Custom imports
from easy_md.utils import fileparser, simulation_util, dcd_image
from easy_md.utils.fileparser import time_tracker

# --------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------
def next_free_state(path: str) -> str:
    """Return a file name with an incremented _N postfix.

    Rules
    -----
    1.  If <path> has no postfix         -> create _0, _1, … until free.
    2.  If <path> ends in _N (integer)   -> bump N until the name is free.
    """
    p      = Path(path)
    parent = p.parent
    parent.mkdir(parents=True, exist_ok=True)      # ensure folder exists

    stem, suffix = p.stem, p.suffix               # 'md_state_1', '.xml'

    # Split at the *last* underscore
    if '_' in stem and stem.rsplit('_', 1)[1].isdigit():
        core, num = stem.rsplit('_', 1)
        num = int(num) + 1                        # start with next integer
    else:
        core, num = stem, 0                       # start fresh at _0

    # Bump until we find a non-existing file
    while True:
        new_path = parent / f'{core}_{num}{suffix}'
        if not new_path.exists():
            return str(new_path)
        num += 1

# --------------------------------------------------------------------------
# Main Simulation Setup
# --------------------------------------------------------------------------
@time_tracker
def main(config=None, starting_state_path=None, starting_checkpoint_path=None, equilibration_only=False):
    """Runs NVT equilibration with monitoring of convergence."""

    # --------------------------------------------------------------------------
    # Setup simulation
    # --------------------------------------------------------------------------
    omm_system, omm_top, off_top = fileparser.load_files(config['path_openmm_system'], config['path_openmm_topology'])
    # simulation_util.print_constraint_info(omm_system, omm_top) # Uncomment to print constraint information

    if config['md_npt']:
        # Set up barostat using flat config structure
        barostat = simulation_util.setup_barostat(
            config['integrator_temperature'],
            config['md_pressure'],
            config['md_barostat_freq'],
            config['md_anisotropic']
        )
        omm_system.addForce(barostat)

    if config['md_harmonic_restraint']:
        force_restraints = simulation_util.setup_force_restraints(reference_structure=config['path_emin_structure'], 
                                                                residue_indices=config['md_restrained_residues'])
        omm_system.addForce(force_restraints)

    # Set up simulation using flat config structure
    simulation = simulation_util.setup_simulation(
        omm_system, 
        omm_top, 
        config['platform_name'],
        {'Precision': config['platform_precision']},
        config['integrator_temperature'],
        config['integrator_friction'],
        config['integrator_timestep']
    )

    # --------------------------------------------------------------------------
    # Load state or checkpoint and setup reporters
    # --------------------------------------------------------------------------
    if starting_state_path is None and starting_checkpoint_path is None:
        starting_state_path = config['path_emin_state']
        print(f"No starting state or checkpoint provided. Using {starting_state_path}")
    
    simulation = simulation_util.load_state_or_checkpoint(
        simulation, 
        temp=config['integrator_temperature'], 
        state_file=starting_state_path,
        checkpoint_file=starting_checkpoint_path
    )

    path_md_image = next_free_state(config['path_md_image'])
    path_md_trajectory = next_free_state(config['path_md_trajectory'])
    path_md_checkpoint = next_free_state(config['path_md_checkpoint'])
    path_md_log = next_free_state(config['path_md_log'])
    

    # Set up reporters
    simulation_util.setup_reporters(
        simulation,
        path_md_log,
        path_md_trajectory,
        path_md_checkpoint,
        config['md_save_interval'],
        config['md_steps']
    )

    # --------------------------------------------------------------------------
    # Run Equilibration
    # --------------------------------------------------------------------------
 
    # Initialize monitoring queues
    temp_window = deque(maxlen=config['monitor_window'])
    energy_window = deque(maxlen=config['monitor_window'])

    if equilibration_only:
        print("\n=== Equilibration ===")
        for step in range(0, config['md_steps'], config['md_save_interval']):
            simulation.step(config['md_save_interval'])

            # Stops equilibration if temperature and energy are within thresholds
            if simulation_util.check_equilibration(
                simulation,
                config['monitor_temp_threshold'],
                config['monitor_energy_threshold'],
                temp_window,
                energy_window,
                config['monitor_window']
            ):
                break
    else:
        print("\n=== Simulation ===")
        simulation.step(config['md_steps'])

    print("Done! Saving state and image")
    simulation.saveState(path_md_image)
    dcd_image.image_molecules(path_md_trajectory, config['path_openff_interchange'], path_md_image)

    print(f"Done! File saved to {next_free_state(config['path_md_state'])}")
    print(f"Trajectory saved to {next_free_state(config['path_md_trajectory'])}")
    print(f"Checkpoint saved to {next_free_state(config['path_md_checkpoint'])}")
    print(f"Log saved to {next_free_state(config['path_md_log'])}")

# --------------------------------------------------------------------------
# Main Execution
# --------------------------------------------------------------------------
if __name__ == '__main__':
    main()