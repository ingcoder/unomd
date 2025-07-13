# Standard library imports
import numpy as np
import os

# Third-party imports
from openmm import(
    LangevinMiddleIntegrator,
    MonteCarloBarostat, 
    MonteCarloAnisotropicBarostat,
    Platform,
    CustomExternalForce,
    openmm,
    unit
)
from openmm.unit import(kilojoule_per_mole, nanometer)
from openmm.app import(Simulation,
                        StateDataReporter, 
                        DCDReporter, 
                        PDBFile,
                        NoCutoff,
                        HBonds,
                        ForceField, 
                        PDBFile
)

# Custom imports
from easy_md.utils import info_logger
import logging

logger = logging.getLogger(__name__)


#------------------------------------------------------------------------------
# Helper Classes
#------------------------------------------------------------------------------

class CheckpointReporter:
    """Reporter for saving simulation checkpoints."""
    def __init__(self, file, reportInterval):
        self._reportInterval = reportInterval
        self._out = open(file, 'wb')

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        simulation.saveCheckpoint(self._out.name)


#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------

def setup_reporters(simulation, log_output, trajectory_output, checkpoint_output, saving_steps, total_steps):
    """Configures and adds state and trajectory reporters to the simulation."""
    # Clear any existing reporters
    simulation.reporters.clear()
    
    reporter_args = {
        'reportInterval': saving_steps,
        'step': True,
        'time': True,
        'potentialEnergy': True,
        'temperature': True,
        'volume': True,
        'separator': '\t',
        'totalSteps': total_steps,
        'speed': True,
        'progress': True
    }
    
    # Add reporters for file and stdout
    for output in [log_output, None]:  # None for stdout
        simulation.reporters.append(StateDataReporter(output, **reporter_args))
    
    simulation.reporters.append(DCDReporter(trajectory_output, saving_steps))

    if checkpoint_output is not None:
        simulation.reporters.append(CheckpointReporter(checkpoint_output, saving_steps))


def setup_simulation(omm_system, omm_top, platform_name, platform_properties, temperature, friction_coef, timestep):
    """Sets up the simulation with specified parameters and returns the configured simulation object."""
    platform = Platform.getPlatformByName(platform_name)
    logger.info(f"Platform being used: {platform.getName()}")

    # Initialize integrator
    integrator = LangevinMiddleIntegrator(
        temperature,
        friction_coef,
        timestep
    )

    # Create simulation
    if platform_name == "CUDA":
        simulation = Simulation(omm_top, omm_system, integrator, platform, platform_properties)
        precision = platform.getPropertyValue(simulation.context, "Precision")
        logger.info(f"Precision being used: {precision}")
    else:
        simulation = Simulation(omm_top, omm_system, integrator, platform)
    
    return simulation

def setup_force_restraints(reference_structure, residue_indices, force_constant=100):
    """Sets up force restraints on the simulation."""
    ref_pdb = PDBFile(reference_structure)
    ref_positions = ref_pdb.getPositions()

    logger.info(f"Adding harmonic positional restraints with force constant {force_constant} kJ/mol/nm^2...")
    # Reduce force constant to avoid instability
    force_k = force_constant * kilojoule_per_mole/nanometer**2   
    restraint_force = CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint_force.addGlobalParameter("k", force_k)
    
    # Add per-particle parameters
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    # Add per-particle reference coordinates
    num_restrained = 0
    chain_restraints = {}

    # Single loop handling both cases
    for idx, atom in enumerate(ref_pdb.topology.atoms()):
        # Combine conditions: atom must be CA and either residue_indices is empty (apply to all) 
        # or the residue id is within the specified range
        if (atom.name == "CA" and 
            (len(residue_indices) == 0 or 
             (residue_indices[0] <= int(atom.residue.id) <= residue_indices[1]))):
            
            x, y, z = ref_positions[idx].value_in_unit(nanometer)
            restraint_force.addParticle(idx, (x, y, z))
            num_restrained += 1
            
            # Update chain restraints using dict.get()
            chain_id = atom.residue.chain.id
            chain_restraints[chain_id] = chain_restraints.get(chain_id, 0) + 1

    logger.info(f"Added positional restraints to {num_restrained} heavy atoms")
    logger.info(f"Restraints per chain: {chain_restraints}")

    return restraint_force

def setup_barostat(temperature, pressure, barostat_frequency, use_anisotropic=False):
    """Configure and return appropriate barostat based on simulation type."""
    logger.info(f"Setting up barostat with pressure {pressure} bar and temperature {temperature} K")
    if use_anisotropic:
        pressure_tuple = (pressure, pressure, pressure)
        scaleXYZ = (True, True, True)
        return MonteCarloAnisotropicBarostat(
            pressure_tuple, temperature,
            scaleXYZ[0], scaleXYZ[1], scaleXYZ[2],
            barostat_frequency
        )
    else:

        return MonteCarloBarostat(pressure, temperature, barostat_frequency)
    
def load_state_or_checkpoint(simulation, temp, state_file=None, checkpoint_file=None):
    """Loads a state or checkpoint file into the simulation to continue equilibration or simulation."""
    try:
        # If load_from_state is explicitly True, or if state_file is provided and exists
        if (state_file and os.path.exists(state_file)):
            if not state_file:
                logger.error("state_file must be provided when load_from_state is True")
                raise ValueError("state_file must be provided when load_from_state is True")
            if not isinstance(state_file, str):
                logger.error("state_file must be a string path")
                raise TypeError("state_file must be a string path")
            
            try:
                simulation.loadState(state_file)
                simulation.currentStep = 0
                simulation.context.setTime(0)
                simulation.context.setVelocitiesToTemperature(temp)
                logger.info(f"Successfully loaded state from {state_file}")
            except Exception as e:
                logger.error(f"Failed to load state file {state_file}: {str(e)}")
                raise RuntimeError(f"Failed to load state file {state_file}: {str(e)}")
                
        elif checkpoint_file:
            if not isinstance(checkpoint_file, str):
                logger.error("checkpoint_file must be a string path")
                raise TypeError("checkpoint_file must be a string path")
                
            try:
                simulation.loadCheckpoint(checkpoint_file)
                logger.info(f"Successfully loaded checkpoint from {checkpoint_file}")
            except Exception as e:
                logger.error(f"Failed to load checkpoint file {checkpoint_file}: {str(e)}")
                raise RuntimeError(f"Failed to load checkpoint file {checkpoint_file}: {str(e)}")
        else:
            logger.error("Either state_file (with load_from_state=True) or checkpoint_file must be provided")
            raise ValueError("Either state_file (with load_from_state=True) or checkpoint_file must be provided")
            
    except FileNotFoundError as e:
        logger.error(f"Error: File not found - {str(e)}")
        raise
    except (ValueError, TypeError) as e:
        logger.error(f"Error: Invalid input - {str(e)}")
        raise
    except RuntimeError as e:
        logger.error(f"Error: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error occurred: {str(e)}")
        raise
        
    return simulation

def print_constraint_info(system, top):
    """Prints constraint information from the system."""
    # Debug: Print topology constraint information
    logger.info(f"Number of constraints: {system.getNumConstraints()}")
    for i in range(min(5, system.getNumConstraints())):  # Print first 5 constraints as example
        constraint = system.getConstraintParameters(i)
        logger.info(f"Constraint {i}: particles {constraint[0]}-{constraint[1]}, distance {constraint[2]}")

    bonds = list(top.bonds())
    logger.info(f"Number of bonds in topology: {len(bonds)}")
    if bonds:
        logger.info("First few bonds:")
        for bond in bonds[:5]:
            logger.info(f"Bond: {bond[0].name}-{bond[1].name}")

def print_system_charge(pdb_filepath, config):
    pdb = PDBFile(pdb_filepath)
    forcefield = ForceField(config['forcefield_protein'], config['forcefield_solvent'])
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
   
    # Iterate over all forces to find the NonbondedForce, which contains charge information
    total_charge = 0.0
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(force.getNumParticles()):
                # For each particle (atom) in the force, extract the charge (the first element of the tuple returned by getParticleParameters)
                charge, _, _ = force.getParticleParameters(i)
                total_charge += charge.value_in_unit(unit.elementary_charge)

    logger.info(f"Total system charge: {total_charge} e")

def check_equilibration(simulation, temp_std_threshold, energy_std_threshold, temp_window, energy_window, window_size):
    """Checks if the system has reached equilibrium based on temperature and energy fluctuations."""
    potential_energy, temperature = get_state_info(simulation)
    # Monitor equilibration
    temp_window.append(temperature)
    energy_window.append(potential_energy)
                                                   
    if len(temp_window) == window_size:
        temp_std = np.std(temp_window)
        energy_std = np.std(energy_window)
        logger.info(f"Temperature std dev: {temp_std:.2f} K")
        logger.info(f"Potential Energy std dev: {energy_std:.2f} kJ/mol")
        if temp_std < temp_std_threshold and energy_std < energy_std_threshold:
            logger.info("\nSystem has reached equilibrium!")
            logger.info(f"Temperature std dev: {temp_std:.2f} K (threshold: {temp_std_threshold} K)")
            logger.info(f"Potential Energy std dev: {energy_std:.2f} kJ/mol (threshold: {energy_std_threshold} kJ/mol)")
            return True
        
    return False

def get_state_info(simulation):
    """Get state information from the simulation."""
    state = simulation.context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    dof = 3 * simulation.system.getNumParticles() - simulation.system.getNumConstraints()
    temperature = (2*state.getKineticEnergy()/(dof*unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)
    return potential_energy, temperature

def transfer_state_with_precision(temp_omm_system, temp_omm_top, checkpoint_file, config):
    """
    If you continue a simulation from a checkpoint that used another precision than the one you want to use.
    First determines the source checkpoint precision by trying different precisions. Then get the state
    from the checkpoint which you can apply to a new simulation object with any precision.
    
    Parameters:
    -----------
    checkpoint_file : str
        Path to the checkpoint file
    omm_top : object
        OpenMM topology object
    system : object
        OpenMM system object
    target_precision : str
        Target precision for the new simulation ('mixed', 'single', or 'double')
    
    Returns:
    --------
    Simulation
        New simulation object with the loaded state and specified precision

    return state

    Example:
        # Create new simulation with target precision
        integrator = LangevinIntegrator(temperature, friction, time_step)
        integrator.setRandomNumberSeed(random_seed)
        platform_properties = {'Precision': target_precision}
        simulation = Simulation(omm_top, system, integrator, platform, platform_properties)
        
        # Transfer state to new simulation
        state = transfer_state_with_precision(....)
        simulation.context.setState(state)
    """
    # Try loading with different precisions to determine the source precision
    precisions = ['double', 'mixed', 'single']
    source_precision = None
    
    for precision in precisions:
        try:
            # We need to create a simulation object with the same precision and platform properties as the checkpoint file. 
            # For example if we ran our first simulation on a CPU with mixed precision and now.
            # Checkpoint does not store precision information we have to use trial and error. Provide the settings for the checkpoint.
            # want to convert it to GPU with double precision. The temporation simulation 
            # has to match the first simulation (CPU, mixed precision)
            temp_integrator = LangevinMiddleIntegrator(config.get('integrator_temperature'), 
                                                       config.get('integrator_friction'), 
                                                       config.get('integrator_timestep'))
            temp_integrator.setRandomNumberSeed(config.get('random_seed'))
            temp_platform = Platform.getPlatformByName('CPU')
            temp_platform_properties = {'Precision': precision}
            temp_simulation = Simulation(temp_omm_top, temp_omm_system, temp_integrator, temp_platform, temp_platform_properties)
            
            # Load the checkpoint with old precision into simulation object
            # Determine simulation 
            temp_simulation.loadCheckpoint(checkpoint_file)
            source_precision = precision
            logger.info(f"Successfully loaded checkpoint with {precision} precision")

            #  Once you have extracted the state from the checkpoint using your 
            # transfer_state_with_precision function, you can apply that state to a 
            # new simulation object with any precision. 
            # The state object contains the numerical data (positions, velocities, etc.) 
            # and is precision-agnostic.
            state = temp_simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            return state
        except Exception as e:
            logger.error(f"Failed to load checkpoint with {precision} precision: {str(e)}")
            continue
    
    if source_precision is None:
        logger.error("Could not determine the precision of the checkpoint file")
        raise ValueError("Could not determine the precision of the checkpoint file")
    