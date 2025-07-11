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

from easy_md.utils.fileparser import load_openff_topology_from_json
import parmed
import os
from pathlib import Path
from typing import Union
import logging
from openff.toolkit import ForceField

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def validate_output_path(file_path: Union[str, Path]) -> Path:
    """Validate that an output file path is writable."""
    path = Path(file_path)
    parent_dir = path.parent
    parent_dir.mkdir(parents=True, exist_ok=True)
    assert os.access(parent_dir, os.W_OK), f"No write permission for directory: {parent_dir}"
    return path


def validate_prmtop(prmtop_path):
    """Validate the AMBER topology file."""
     
    # Load the AMBER files into a ParmEd structure for validation
    #prmtop_path needs to be string.
    amber_structure = parmed.load_file(prmtop_path)

    # # Validate input files exist and are readable
    # prmtop_path = Path(prmtop_path)
    # # inpcrd_path = Path(inpcrd_path)
    # assert prmtop_path.exists(), f"Topology file not found: {prmtop_path}"
    # # assert inpcrd_path.exists(), f"Coordinate file not found: {inpcrd_path}"
    # assert os.access(prmtop_path, os.R_OK), f"No read permission for topology file: {prmtop_path}"
    # # assert os.access(inpcrd_path, os.R_OK), f"No read permission for coordinate file: {inpcrd_path}"
   

    # Check unique residues 
    unique_residues = set(residue.name for residue in amber_structure.residues)
    logger.info(f"Residues in AMBER structure: {unique_residues}")

    # Check for water and ions and ligand
    water_residues = [res for res in amber_structure.residues if res.name in ['WAT', 'HOH']]
    sodium_ions = [res for res in amber_structure.residues if res.name == 'NA']
    chloride_ions = [res for res in amber_structure.residues if res.name == 'CL']
    ligand = [res for res in amber_structure.residues if res.name in ['UNK', 'LIG']]
    logger.info(f"Number of water molecules: {len(water_residues)}")
    logger.info(f"Number of sodium ions: {len(sodium_ions)}")
    logger.info(f"Number of chloride ions: {len(chloride_ions)}")
    logger.info(f"Number of ligands: {len(ligand)}")

    # Check for multiple molecules
    logger.info("\nChecking for multiple molecules:")
    molecules = []
    current_molecule = set()
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
    logger.info(f"Found {len(molecules)} separate molecules")

    # Check bonds with missing type information
    problematic_bonds = []
    for bond in amber_structure.bonds:
        if bond.type is None:
            atom1_info = f"{bond.atom1.residue.name} {bond.atom1.residue.idx+1} {bond.atom1.name}"
            atom2_info = f"{bond.atom2.residue.name} {bond.atom2.residue.idx+1} {bond.atom2.name}"
            problematic_bonds.append(f"Bond between {atom1_info} and {atom2_info}")
    
    if problematic_bonds:
        logger.warning("⚠️ Warning: Found bonds with missing type information:")
        for bond in problematic_bonds:
            logger.warning(f"  - {bond}")
        logger.warning("\nThis might indicate:")
        logger.warning("1. Missing force field parameters")
        logger.warning("2. Incompatible force field types")
        logger.warning("3. Issues with ligand parameters")
        raise ValueError("Cannot create Amber topology due to missing bond type information")
    
    logger.info("✅ Successfully validated AMBER topology file")


#------------------------------------------------------------------------------
# Conversion
#------------------------------------------------------------------------------
def omm_to_prmtop(off_top, output_prmtop_path, output_inpcrd_path):
    """Convert OpenFF topology to AMBER topology."""
    # Validate input file
    off_top = Path(off_top)
    assert off_top.exists(), f"Input file not found: {off_top}"
    assert os.access(off_top, os.R_OK), f"No read permission for file: {off_top}"
    
    # Validate output paths
    output_prmtop_path = validate_output_path(output_prmtop_path)
    output_inpcrd_path = validate_output_path(output_inpcrd_path)

    off_top = load_openff_topology_from_json(off_top)
    sage_ff14sb = ForceField("openff-2.0.0.offxml", "ff14sb_off_impropers_0.0.0.3.offxml")
    interchange = sage_ff14sb.create_interchange(off_top)
    omm_system = interchange.to_openmm(add_constrained_forces=True)
    omm_top = interchange.to_openmm_topology()
    parmed_structure = parmed.openmm.load_topology(omm_top, omm_system)
    parmed_structure.save(output_prmtop_path, overwrite=True)


def omm_to_prmtop2(off_top, output_prmtop_path):
    # Convert to ParmEd Structure
    # omm_system = load_openmm_system_from_xml(omm_system)
    # omm_top = load_openmm_topology_from_pickle(omm_top)
    off_top = load_openff_topology_from_json(off_top)

    sage_ff14sb = ForceField("openff-2.0.0.offxml", "ff14sb_off_impropers_0.0.3.offxml")
    interchange = sage_ff14sb.create_interchange(off_top)

    omm_system = interchange.to_openmm(add_constrained_forces=True)
    omm_top = interchange.to_openmm_topology()

    # Create ParmEd structure
    parmed_structure = parmed.openmm.load_topology(omm_top, omm_system)
    parmed_structure.save(output_prmtop_path, overwrite=True)



#------------------------------------------------------------------------------
# Main Function
#------------------------------------------------------------------------------
def main():
    try:
        # Define paths
        off_topology_json_path = '/Users/ingrid/Projects/EasyMD/easy-md/example/output/openff_topology.json'
        output_prmtop_path = '/Users/ingrid/Projects/EasyMD/easy-md/example/output/4W52_solvated_complex.prmtop'
        output_inpcrd_path = '/Users/ingrid/Projects/EasyMD/easy-md/example/output/4W52_solvated_complex.inpcrd'

        # Convert OpenFF topology to AMBER topology
        omm_to_prmtop2(off_topology_json_path, output_prmtop_path)
        validate_prmtop(output_prmtop_path)
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