# from unomd.utils.config import create_config
# from unomd.main import run_receptor_preperation

# config = create_config(
#     protein_file="_TRIMED.pdb"
# )

# run_receptor_preperation.prepare_receptor(config)


from unomd.main.quickrun import quickrun

# Run a simple protein-ligand simulation
quickrun(
    protein_file="TRIMED_final.pdb",
    nsteps=1000
)