from easy_md.main import run_solvation, run_forcefield_parameterization, run_energy_minimization, run_simulation
from easy_md.utils.config import create_config

def quickrun(protein_file, ligand_file=None, nsteps=1000):
    config = create_config(
        protein_file=protein_file,
        ligand_file=ligand_file,
        md_steps=nsteps
    )
    run_solvation.add_water(config=config)
    run_forcefield_parameterization.main(config)
    run_energy_minimization.main(config)
    run_simulation.main(config)