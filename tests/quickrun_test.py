from unomd.main.quickrun import quickrun

# Run a simple protein-ligand simulation
quickrun(
    protein_file="/Users/ingrid/Projects/UnoMD/unomd/example/4W52.pdb",
    ligand_file="/Users/ingrid/Projects/UnoMD/unomd/example/4w52_C_EPE.sdf",
    nsteps=100,
    export_prmtop=True,
    mmpbsa=True
)