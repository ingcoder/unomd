"""
This script loads a PDB file, calculates the bounding box dimensions, and solvates the system using OpenMM.

### Overview:
1. **Load PDB file**: Reads the protein structure from a PDB file.
2. **Compute bounding box**: Determines the dimensions of the protein.
3. **Define simulation box**: Adds a buffer and sets up box vectors.
4. **Add solvent**: Uses OpenMM's `Modeller` to solvate the system.
5. **Save solvated structure**: Outputs the solvated system to a PDB file.

### Parameters:
- `PDB_FILE`: Path to the input PDB file.
- `SOLVATED_FILE`: Path to save the solvated PDB file.
- `BUFFER`: Extra space added to the bounding box for solvation (default: 2.5 nm).
- `IONIC_STRENGTH`: The ionic concentration for solvation (default: 0.15 M).
- `FORCEFIELD_FILES`: Force field files used for modeling the system.

### Usage:
Run the script:
```bash
python script.py
"""

# Standard library imports
import numpy as np
import subprocess

# Third-party imports
from openmm.app import PDBFile, Modeller, ForceField
from openmm import Vec3
from openmm.unit import nanometer as nm, molar
from pdbfixer import PDBFixer

# Custom imports & logging
from unomd.utils.fileparser import time_tracker
from unomd.utils import info_logger
import logging


logger = logging.getLogger(__name__)

def run_cmd(command):
    """Execute a shell command and handle errors."""
    try:
        logger.info(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True)
        logger.info("✅ Command executed successfully")
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Command failed with return code {e.returncode}: {command}")
        raise
    except FileNotFoundError:
        logger.error(f"❌ Command not found: {command}")
        raise

def convert_pdb_to_pqr(input_file, output_file):
    run_cmd(f"pdb2pqr --ff AMBER --with-ph=7.4 {input_file} {output_file}")

def convert_pqr_to_pdb(input_file, output_file):
    run_cmd(f"obabel {input_file} -O {output_file}")


@time_tracker
def prepare_receptor(config):
    logger.info("========================================================")
    logger.info(" Receptor Preparation")
    logger.info("========================================================")

    logger.info("Loading PDB file & adding missing hydrogens...")
    fixer = PDBFixer(filename=config.get("path_protein"))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    n_missing_heavy = sum(len(v) for v in fixer.missingAtoms.values())

    if n_missing_heavy > 0:
        logger.info(f"Found {n_missing_heavy} missing heavy atoms - adding them now...")
        fixer.addMissingAtoms()
        # logger.info("Adding missing hydrogens...")
        # fixer.addMissingHydrogens(pH=config.get("solv_pH"))
    else:
        logger.info("No missing heavy atoms found")
        #logger.info("Adding missing hydrogens...")
        # fixer.addMissingHydrogens(pH=config.get("solv_pH"))

    # Save the processed fixer output to a temporary PDB file
    temp_processed_pdb = config.get("path_protein").replace(".pdb", "_processed.pdb")
    logger.info(f"Saving processed structure from fixer to: {temp_processed_pdb}")
    with open(temp_processed_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # Now use the processed PDB for the conversion chain
    logger.info("Converting processed PDB to PQR...")
    temp_protein_pqr = config.get("path_protein").replace(".pdb", ".pqr")
    convert_pdb_to_pqr(temp_processed_pdb, temp_protein_pqr)

    logger.info("Converting PQR back to final PDB...")
    final_pdb = config.get("path_protein").replace(".pdb", "_final.pdb")
    convert_pqr_to_pdb(temp_protein_pqr, final_pdb)
    
    logger.info(f"Final processed PDB saved to: {final_pdb}")
    logger.info("✅ Receptor preparation completed successfully")

if __name__ == "__main__":
    prepare_receptor()
