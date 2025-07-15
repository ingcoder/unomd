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
import yaml
import glob

# Third-party imports
import parmed
from openff.toolkit import ForceField

# Custom imports
from easy_md.utils.fileparser import load_openff_topology_from_json, time_tracker
from easy_md.utils import info_logger
import logging

logger = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

WATER_RESIDUES = ["WAT", "HOH"]
LIGAND_RESIDUES = ["UNK", "LIG"]
ION_RESIDUES = {"NA": "sodium", "CL": "chloride"}


# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------


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

    logger.info(f"Residue counts: {residue_counts}")

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
    logger.info(f"Number of seperate molecules: {len(molecules)}")

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

    logger.info(f"Problematic bonds: {problematic_bonds}")
    return problematic_bonds


def validate_prmtop(config) -> None:
    logger.info("Validating AMBER topology file...")
    amber_structure = parmed.load_file(config["path_amber_solvated"])

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


def get_last_final_dcd(config):
    """Getting the last final dcd file."""
    try:
        trajectory_dir = os.path.dirname(config["path_md_image"])
        dcd_files = glob.glob(os.path.join(trajectory_dir, "final_md*.dcd"))
        dcd_files.sort(key=os.path.getmtime)
        logger.info(f"Last final dcd file: {dcd_files[-1]}")
        return dcd_files[-1]
    except Exception as e:
        logger.error(f"Error getting last final dcd file: {e}")
        raise 
    
def clean_mmpbsa_files():
    """Cleaning the mmpbsa files."""
    script_dir = Path(__file__).resolve().parent
    for file in script_dir.glob("_MMPBSA*"):
        try:
            file.unlink()
            logger.info(f"Deleted: {file}")
        except Exception as e:
            logger.error(f"Could not delete {file}: {e}")


# ------------------------------------------------------------------------------
# Conversion Functions
# ------------------------------------------------------------------------------


def offtopology_to_prmtop(config) -> None:
    logger.info(f"Converting OpenFF topology to AMBER topology...")

    # Load and convert topology
    off_topology = load_openff_topology_from_json(config["path_openff_topology"])

    found = False
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

    # sage_ff14sb = ForceField("openff-2.0.0.offxml", "ff14sb_off_impropers_0.0.3.offxml")
    sage_ff14sb = ForceField(
        config["ff_small_molecule_openff"], config["ff_protein_openff"]
    )
    interchange = sage_ff14sb.create_interchange(off_topology)

    # Convert to OpenMM
    omm_system = interchange.to_openmm(add_constrained_forces=True)
    omm_topology = interchange.to_openmm_topology()

    # Create ParmEd structure and save
    parmed_structure = parmed.openmm.load_topology(omm_topology, omm_system)
    parmed_structure.save(str(config["path_amber_solvated"]), overwrite=True)

    logger.info(
        f"✅ Successfully converted OpenFF topology to AMBER topology: {config['path_amber_solvated']}"
    )

@time_tracker
def run_ante_mmpbsa(config) -> None:
    """Preparing the amber input topology files for the MMPBSA calculation."""

    command = [
        "ante-MMPBSA.py",
        "-p",
        config["path_amber_solvated"],
        "-c",
        config["path_amber_complex"],
        "-r",
        config["path_amber_receptor"],
        "-l",
        config["path_amber_ligand"],
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

@time_tracker
def run_mmpbsa(config) -> None:
    """Running the MMPBSA calculation."""

    last_final_dcd = get_last_final_dcd(config)

    command = [
        "MMPBSA.py",
        "-O",
        "-i",
        config["path_mmpbsa_in"],
        "-o",
        config["path_mmpbsa_results"],
        "-sp",
        config["path_amber_solvated"],
        "-cp",
        config["path_amber_complex"],
        "-rp",
        config["path_amber_receptor"],
        "-lp",
        config["path_amber_ligand"],
        "-y",
        last_final_dcd,
    ]

    logger.info("Running mmpbsa.py with the following command:")
    logger.info(" ".join(command))

    try:
        subprocess.run(command, check=True)
        logger.info("✅ mmpbsa.py ran successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ mmpbsa.py failed with return code {e.returncode}")
    except FileNotFoundError:
        logger.error("❌ mmpbsa.py not found. Is it installed and in your PATH?")


def create_mmpbsa_input_file(config):
    """Creating the mmpbsa input file."""
    try:
        with open(config["path_mmpbsa_in"], "w") as f:
            f.write("&general\n")
            f.write("interval   = 10,\n")
            f.write("verbose    = 2,\n")
            f.write("endframe=100,\n")
            f.write("keep_files=2,\n")
            f.write("strip_mask= :HOH:CL:NA,\n")
            f.write("\n")
            f.write("/\n")
            f.write("&gb\n")
            f.write("igb=2, saltcon=0.150,\n")
            f.write("/\n")
            f.write("&pb\n")
            f.write("istrng=0.150,inp=2, radiopt=0, prbrad=1.4,\n")
            f.write("/\n")
        logger.info(f"✅ Successfully created mmpbsa input file: {config['path_mmpbsa_in']}")
    except Exception as e:
        logger.error(f"Error creating mmpbsa input file: {e}")
        raise 


def export_prmtop(config):
    # Convert OpenFF topology to AMBER topology
    offtopology_to_prmtop(config)
    validate_prmtop(config)

# ------------------------------------------------------------------------------
# Main Function
# ------------------------------------------------------------------------------


def main(config_filepath) -> None:
    logger.info("========================================================")
    logger.info("🛠️ MM/PBSA Calculation")
    logger.info("========================================================")

    if config_filepath is not None:
        with open(config_filepath, "r") as f:
            config = yaml.safe_load(f)
    else:
        logger.error("No config file provided. Please provide a config file.")
        raise ValueError("No config file provided. Please provide a config file.")

    create_mmpbsa_input_file(config)

    try:
        export_prmtop(config)
        run_ante_mmpbsa(config)
        run_mmpbsa(config)

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
    config_filepath = (
        "/Users/ingrid/Projects/EasyMD/easy-md/example/config/simulation_config.yaml"
    )
    clean_mmpbsa_files()
    # main(config_filepath)
