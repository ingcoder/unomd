"""
AMBER file utilities for molecular dynamics simulations.

This module provides functions for converting OpenFF topologies to AMBER format,
validating AMBER topology files, and handling residue counting and molecule
detection in AMBER structures.
"""

from typing import List, Set
import parmed
from openff.toolkit import ForceField
from unomd.utils.fileparser import load_openff_topology_from_json
from unomd.utils import info_logger
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


def export_prmtop(config):
    """Export the AMBER topology to a PRM topology file."""
    # Convert OpenFF topology to AMBER topology
    offtopology_to_prmtop(config)
    validate_prmtop(config)
