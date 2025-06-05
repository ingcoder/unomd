from easy_md.main import run_solvation, run_forcefield_parameterization, run_energy_minimization, run_simulation
from easy_md.utils.config import create_config

def quickrun(protein_file, ligand_file=None, nsteps=1000, **kwargs):
    """Run a quick molecular dynamics simulation with customizable parameters.
    
    Args:
        protein_file (str): Path to protein PDB file
        ligand_file (str, optional): Path to ligand file. Defaults to None.
        nsteps (int, optional): Number of MD steps. Defaults to 1000.
        **kwargs: Additional configuration parameters. Can include any parameter supported by create_config:
            
            Platform Settings:
                platform_name (str): Platform for computation. Default: "CPU". Options: ['CPU', 'CUDA']
                platform_precision (str): Precision mode. Default: "mixed". Options: ['single', 'mixed', 'double']
            
            Integrator Settings:
                integrator_temperature (float): Temperature in Kelvin. Default: 300.0
                integrator_friction (float): Friction coefficient in ps^-1. Default: 1.0
                integrator_timestep (float): Time step in ps. Default: 0.002
            
            Solvation Settings:
                solv_box_buffer (float): Buffer size in angstroms. Default: 2.5
                solv_ionic_strength (float): Ionic strength in molar. Default: 0.15
                solv_positive_ion (str): Type of positive ion. Default: "Na+"
                solv_negative_ion (str): Type of negative ion. Default: "Cl-"
                solv_model (str): Water model. Default: "tip3p"
                solv_pH (float): pH of the solvent. Default: 7.0
            
            MD Simulation Settings:
                md_save_interval (int): Save interval for trajectory. Default: 10
                md_pressure (float): Pressure in atmospheres. Default: 1.0
                md_anisotropic (bool): Use anisotropic pressure. Default: False
                md_barostat_freq (int): Barostat frequency. Default: 25
                md_harmonic_restraint (bool): Use harmonic restraints. Default: True
                md_load_state (bool): Load previous state if available. Default: True
                md_restrained_residues (list): List of residues to restrain. Default: []
                md_npt (bool): Use NPT ensemble. Default: False
            
            And many more - see create_config documentation for full list.
    """
    # Create configuration with all provided parameters
    config = create_config(
        protein_file=protein_file,
        ligand_file=ligand_file,
        md_steps=nsteps,
        **kwargs
    )
    
    # Run the simulation pipeline
    run_solvation.add_water(config=config)
    run_forcefield_parameterization.main(config)
    run_energy_minimization.main(config)
    run_simulation.main(config)