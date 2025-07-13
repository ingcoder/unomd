"""Script for converting OpenMM/OpenFF molecular systems to AMBER format for MM/PBSA calculations.

This script provides functionality to:
1. Convert OpenMM/OpenFF molecular systems to AMBER format (prmtop/inpcrd)
2. Validate the resulting AMBER topology files
3. Handle force field parameters and bond types
4. Support both protein-ligand complexes and standalone molecules

The script includes robust error checking for:
- Bonds with missing type information
- Ligand parameter issues
- Counts number of molecules in the system

The conversion process uses OpenFF's Interchange to handle force field parameters
and ParmEd for the final AMBER format conversion. It supports the ff14SB force field
with OpenFF impropers for proteins and SMIRNOFF for ligands.

Note: The script is specifically designed to handle systems that will be used
for MM/PBSA calculations, ensuring proper handling of constrained bonds and
force field parameters.
"""

# Standard library imports
import logging
import os
from pathlib import Path
from typing import List, Set, Union
import subprocess

# Third-party imports
import parmed
from openff.toolkit import ForceField

# Custom imports
from easy_md.utils.fileparser import load_openff_topology_from_json
from easy_md.utils import info_logger
import logging

logger = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

FF_SMALL_MOLECULE = ("openff-2.0.0.offxml",)
FF_PROTEIN = "ff14sb_off_impropers_0.0.0.3.offxml"
WATER_RESIDUES = ["WAT", "HOH"]
LIGAND_RESIDUES = ["UNK", "LIG"]
ION_RESIDUES = {"NA": "sodium", "CL": "chloride"}
DEFAULT_INPUT_PATH = (
    "/Users/ingrid/Projects/EasyMD/easy-md/example/output/openff_topology.json"
)
DEFAULT_OUTPUT_PATH = (
    "/Users/ingrid/Projects/EasyMD/easy-md/example/output/4W52_solvated_complex.prmtop"
)


# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------


def validate_output_path(file_path: Union[str, Path]) -> Path:
    path = Path(file_path)
    parent_dir = path.parent
    parent_dir.mkdir(parents=True, exist_ok=True)

    if not os.access(parent_dir, os.W_OK):
        logger.error(f"No write permission for directory: {parent_dir}")
        raise AssertionError(f"No write permission for directory: {parent_dir}")

    return path


def _count_residues_by_type(amber_structure: parmed.Structure) -> dict:
    residue_counts = {}

    # Count water molecules
    water_residues = [
        res for res in amber_structure.residues if res.name in WATER_RESIDUES
    ]
    residue_counts["water"] = len(water_residues)

    # Count ions
    for ion_code, ion_name in ION_RESIDUES.items():
        ions = [res for res in amber_structure.residues if res.name == ion_code]
        residue_counts[ion_name] = len(ions)

    # Count ligands
    ligand_residues = [
        res for res in amber_structure.residues if res.name in LIGAND_RESIDUES
    ]
    residue_counts["ligand"] = len(ligand_residues)

    return residue_counts


def _find_molecules(amber_structure: parmed.Structure) -> List[Set]:
    """Find separate molecules in the AMBER structure."""
    molecules = []
    remaining_atoms = set(amber_structure.atoms)

    while remaining_atoms:
        start_atom = remaining_atoms.pop()
        current_molecule = {start_atom}
        to_process = {start_atom}

        while to_process:
            atom = to_process.pop()
            for bonded_atom in atom.bond_partners:
                if bonded_atom in remaining_atoms:
                    current_molecule.add(bonded_atom)
                    to_process.add(bonded_atom)
                    remaining_atoms.remove(bonded_atom)

        molecules.append(current_molecule)

    return molecules


def _check_problematic_bonds(amber_structure: parmed.Structure) -> List[str]:
    """Check for bonds with missing type information."""
    problematic_bonds = []

    for bond in amber_structure.bonds:
        if bond.type is None:
            atom1_info = (
                f"{bond.atom1.residue.name} "
                f"{bond.atom1.residue.idx + 1} "
                f"{bond.atom1.name}"
            )
            atom2_info = (
                f"{bond.atom2.residue.name} "
                f"{bond.atom2.residue.idx + 1} "
                f"{bond.atom2.name}"
            )
            problematic_bonds.append(f"Bond between {atom1_info} and {atom2_info}")

    return problematic_bonds


def validate_prmtop(prmtop_path: Union[str, Path]) -> None:
    logger.info("Validating AMBER topology file...")
    amber_structure = parmed.load_file(str(prmtop_path))

    # Check unique residues
    unique_residues = set(residue.name for residue in amber_structure.residues)
    logger.info(f"Residues in AMBER structure: {unique_residues}")

    # Count different types of residues
    residue_counts = _count_residues_by_type(amber_structure)
    logger.info(f"Number of water molecules: {residue_counts['water']}")
    logger.info(f"Number of sodium ions: {residue_counts['sodium']}")
    logger.info(f"Number of chloride ions: {residue_counts['chloride']}")
    logger.info(f"Number of ligands: {residue_counts['ligand']}")

    # Check for multiple molecules
    logger.info("Checking for multiple molecules:")
    molecules = _find_molecules(amber_structure)
    logger.info(f"Found {len(molecules)} separate molecules")

    # Check bonds with missing type information
    problematic_bonds = _check_problematic_bonds(amber_structure)

    if problematic_bonds:
        logger.warning("⚠️ Warning: Found bonds with missing type information:")
        for bond in problematic_bonds:
            logger.warning(f"  - {bond}")
        logger.warning("This might indicate:")
        logger.warning("1. Missing force field parameters")
        logger.warning("2. Incompatible force field types")
        logger.warning("3. Issues with ligand parameters")
        raise ValueError(
            "Cannot create Amber topology due to missing bond type information"
        )

    logger.info("✅ Successfully validated AMBER topology file")


# ------------------------------------------------------------------------------
# Conversion Functions
# ------------------------------------------------------------------------------


def offtopology_to_prmtop(
    off_top: Union[str, Path], output_prmtop_path: Union[str, Path]
) -> None:
    logger.info(f"Converting OpenFF topology to AMBER topology...")

    # Load and convert topology
    off_topology = load_openff_topology_from_json(off_top)

    found = False
    unique_residues = set()
    for molecule in off_topology.molecules:
        for atom in molecule.atoms:
            atom_count = molecule.n_atoms
            if atom_count < 100:  # heuristic: ligands often <100 atoms
                found = True
                break

    if found:
        logger.info(f"✅ Found ligand in topology (molecule under 100 atoms)")
    else:
        logger.error("❌ ligand not found in topology. Halting conversion.")
        raise ValueError("ligand not found in the OpenFF topology.")

    sage_ff14sb = ForceField("openff-2.0.0.offxml", "ff14sb_off_impropers_0.0.3.offxml")
    # sage_ff14sb = ForceField(FF_SMALL_MOLECULE, FF_PROTEIN)
    interchange = sage_ff14sb.create_interchange(off_topology)

    # Convert to OpenMM
    omm_system = interchange.to_openmm(add_constrained_forces=True)
    omm_topology = interchange.to_openmm_topology()

    # Create ParmEd structure and save
    parmed_structure = parmed.openmm.load_topology(omm_topology, omm_system)
    parmed_structure.save(str(output_prmtop_path), overwrite=True)

    logger.info(
        f"✅ Successfully converted OpenFF topology to AMBER topology: {output_prmtop_path}"
    )


def run_ante_mmpbsa(prmtop_path: Union[str, Path]) -> None:
    command = [
        "ante-MMPBSA.py",
        "-p",
        prmtop_path,
        "-c",
        "/Users/ingrid/Projects/EasyMD/easy-md/example/output/system_complex.prmtop",
        "-r",
        "/Users/ingrid/Projects/EasyMD/easy-md/example/output/receptor.prmtop",
        "-l",
        "/Users/ingrid/Projects/EasyMD/easy-md/example/output/ligand.prmtop",
        "-s",
        ":HOH:CL:NA",
        "-n",
        ":UNK",
        "--radii",
        "mbondi2",
    ]

    logger.info("Running ante-MMPBSA.py with the following command:")
    logger.info(" ".join(command))

    try:
        subprocess.run(command, check=True)
        logger.info("✅ ante-MMPBSA.py ran successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ ante-MMPBSA.py failed with return code {e.returncode}")
    except FileNotFoundError:
        logger.error("❌ ante-MMPBSA.py not found. Is it installed and in your PATH?")


# ------------------------------------------------------------------------------
# Main Function
# ------------------------------------------------------------------------------


def main() -> None:
    logger.info("========================================================")
    logger.info("🛠️ MM/PBSA Calculation")
    logger.info("========================================================")

    try:
        # Define paths
        off_topology_json_path = DEFAULT_INPUT_PATH
        output_prmtop_path = DEFAULT_OUTPUT_PATH

        # Convert OpenFF topology to AMBER topology
        offtopology_to_prmtop(off_topology_json_path, output_prmtop_path)
        validate_prmtop(output_prmtop_path)
        run_ante_mmpbsa(output_prmtop_path)

    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        raise
    except PermissionError as e:
        logger.error(f"Permission error: {e}")
        raise
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


if __name__ == "__main__":
    main()
