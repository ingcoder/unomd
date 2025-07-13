"""
This script performs energy minimization and equilibration of a molecular system using OpenMM.
It loads a pre-prepared system and topology, performs energy minimization, gradual heating,
equilibration, and saves the final minimized structure as pdb and xml openMM state files.
"""

# Standard library imports
import numpy as np
from openff.toolkit import Topology

# Third-party imports
import openmm
from openmm.app import PDBFile
from openmm.unit import kelvin

# Custom imports
from easy_md.utils import fileparser
from easy_md.utils import simulation_util
from easy_md.utils.fileparser import time_tracker
from easy_md.utils import info_logger
import logging

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------

def energy_force_post_simulation(simulation):
    """Calculates and prints the potential energy and top 10 force magnitudes after simulation."""
    potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    # Calculate the magnitude of the force vectors and find the maximum
    forces = simulation.context.getState(
        getForces=True
    ).getForces()  # This is a list of Vec3 objects
    max_force_magnitude = [(f.x**2 + f.y**2 + f.z**2) ** 0.5 for f in forces]
    top_10_force_magnitudes = np.sort(max_force_magnitude)[-10:]
    logger.info(f"Potential Energy after Minimization: {potential_energy}")
    logger.info(
        f"Maximum Force Magnitude after Minimization: {top_10_force_magnitudes}"
    )


def describe_state(state: openmm.State, name: str = "State"):
    """Prints the potential energy and maximum force for a given state."""
    max_force = max(np.sqrt(v.x**2 + v.y**2 + v.z**2) for v in state.getForces())
    logger.info(
        f"{name} has energy {round(state.getPotentialEnergy()._value, 2)} kJ/mol "
        f"with maximum force {round(max_force, 2)} kJ/(mol nm)"
    )


def save_min_structure(simulation, emin_pdb_output, emin_xml_output):
    """Saves the minimized structure as PDB and XML files."""
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(emin_pdb_output, "w"))
    simulation.saveState(emin_xml_output)


# --------------------------------------------------------------------------
# Main Simulation Setup
# --------------------------------------------------------------------------

@time_tracker
def main(config):
    logger.info("========================================================")
    logger.info("⬇️ Energy Minimization")
    logger.info("========================================================")

    # Load OpenFF topology of solvated system
    off_top = Topology.from_json(open(config["path_openff_topology"]).read())
    omm_system, omm_top, _ = fileparser.load_files(
        config["path_openmm_system"], config["path_openmm_topology"]
    )

    # Set up simulation using flat config structure
    emin_simulation = simulation_util.setup_simulation(
        omm_system,
        omm_top,
        config["platform_name"],
        {"Precision": config["platform_precision"]},
        config["integrator_temperature"],
        config["integrator_friction"],
        config["integrator_timestep"],
    )

    emin_simulation.context.setPositions(off_top.get_positions().to_openmm())

    # Initial minimization before heating
    logger.info("Starting Energy Minimization...")
    emin_simulation.minimizeEnergy()
    initial_state = emin_simulation.context.getState(getEnergy=True, getForces=True)
    describe_state(initial_state, "Initial minimized state")

    # Gradual heating process to prevent energy spikes and to let the system naturally adapt to new temperature.
    logger.info("Gradual Heating Process...")
    for temp in range(0, config["emin_target_temp"] + 1, config["emin_heating_step"]):
        emin_simulation.context.setVelocitiesToTemperature(temp * kelvin)
        emin_simulation.step(config["emin_heating_interval"])
        current_state = emin_simulation.context.getState(getEnergy=True, getForces=True)
        describe_state(current_state, f"Heating at {temp}K")

    # Final equilibration
    logger.info(f"Final Equilibration at {config['emin_target_temp']}K...")
    emin_simulation.step(config["emin_steps"])

    # Final minimization. Adjusts atomic positions to reduce energy without considering time evolution.
    logger.info("Final Energy Minimization...")
    emin_simulation.minimizeEnergy()
    final_state = emin_simulation.context.getState(getEnergy=True, getForces=True)
    describe_state(final_state, "Final minimized state")

    save_min_structure(
        emin_simulation, config["path_emin_structure"], config["path_emin_state"]
    )

    logger.info(
        f"✅ Energy Minimization and Equilibration Complete! File saved to: {config['path_emin_structure']}"
    )


# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
