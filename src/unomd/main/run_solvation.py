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


@time_tracker
def add_water(config):
    logger.info("========================================================")
    logger.info("💧 Solvation")
    logger.info("========================================================")

    logger.info("Loading PDB file & adding missing hydrogens...")
    fixer = PDBFixer(filename=config.get("path_protein"))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    n_missing_heavy = sum(len(v) for v in fixer.missingAtoms.values())

    if n_missing_heavy > 0:
        logger.info(f"Found {n_missing_heavy} missing heavy atoms - adding them now...")
        fixer.addMissingAtoms()
        logger.info("Adding missing hydrogens...")
        fixer.addMissingHydrogens(pH=config.get("solv_pH"))
    else:
        logger.info("No missing heavy atoms found")
        logger.info("Adding missing hydrogens...")
        fixer.addMissingHydrogens(pH=config.get("solv_pH"))

    # Create Modeller instance from fixed structure
    modeller = Modeller(fixer.topology, fixer.positions)

    # Extract positions and convert to numpy array
    positions = np.array(
        [[atom.x, atom.y, atom.z] for atom in modeller.positions.value_in_unit(nm)]
    )

    # Calculate the min and max along each axis
    min_coords = np.min(positions, axis=0)
    max_coords = np.max(positions, axis=0)
    box_dimensions = max_coords - min_coords

    # Get forcefield from config
    forcefield = ForceField(config.get("ff_protein"), config.get("ff_water"))

    # Define box dimensions
    x_dimension = box_dimensions[0] + config.get("solv_box_buffer")
    y_dimension = box_dimensions[1] + config.get("solv_box_buffer")
    z_dimension = box_dimensions[2] + config.get("solv_box_buffer")

    logger.info("Water Box Dimensions (nanometers):")
    logger.info(f"Width (X-axis): {x_dimension}")
    logger.info(f"Height (Y-axis): {y_dimension}")
    logger.info(f"Depth (Z-axis): {z_dimension}")

    box_vecs = (
        Vec3(x_dimension, 0, 0) * nm,
        Vec3(0, y_dimension, 0) * nm,
        Vec3(0, 0, z_dimension) * nm,
    )

    logger.info("Adding water, Na+, Cl- ions...")
    try:
        modeller.addSolvent(
            forcefield,
            boxVectors=box_vecs,
            ionicStrength=config.get("solv_ionic_strength") * molar,
            positiveIon=config.get("solv_positive_ion"),
            negativeIon=config.get("solv_negative_ion"),
            model=config.get("solv_model", "tip3p"),
        )
    except Exception as e:
        logger.error(f"❌ Error adding solvent: {e}")
        raise

    with open(config.get("path_protein_solvated"), "w") as file:
        PDBFile.writeFile(modeller.topology, modeller.positions, file)
    logger.info(
        f"✅ Saved solvated structure to: {config.get('path_protein_solvated')}"
    )


if __name__ == "__main__":
    add_water()
