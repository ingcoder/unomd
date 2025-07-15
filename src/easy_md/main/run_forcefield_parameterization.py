"""
This module handles the creation and parameterization of molecular systems for OpenMM simulations.
It performs the following key tasks:
1. Loads and processes protein and ligand structures combined in a single complex
2. Creates OpenFF topology for the protein-ligand complex
3. Parameterizes the complex system using OpenFF forcefields
4. Converts and validates OpenMM system and topology and saves PDB from OpenFF Topology
"""

# Standard library imports
from typing import Iterable
import pickle
import warnings
import numpy as np

# Third-party imports
from openmm import unit, XmlSerializer
from openmm.app import PDBFile
from openff.toolkit import Molecule, Topology, ForceField
from openff.units import Quantity, unit

# Costum imports
from easy_md.utils.fileparser import time_tracker
from easy_md.utils import ligand_util
from easy_md.utils import info_logger
import logging

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------


def check_for_large_molecules(topology, atom_count_threshold=100):
    """Sanity check if there are any large molecules in the topology,
    which might indicate the presence of proteins or polymers."""
    for molecule in topology.molecules:
        if len(molecule.atoms) > atom_count_threshold:
            logger.info(
                f"Found a large molecule with {len(molecule.atoms)} atoms, which might be a protein or polymer."
            )
            return True
    logger.error("No large molecules found that could indicate a protein or polymer.")
    return False


def load_files(pdb_solv_file, sdf_file):
    try:
        openff_protein_top = Topology.from_pdb(pdb_solv_file)
    except Exception as e:
        logger.error(f"Error loading protein topology: {e}")

    if sdf_file == "":
        openff_ligand_mol = None
    else:
        openff_ligand_mol = Molecule.from_file(sdf_file)

    return openff_protein_top, openff_ligand_mol


def validate_system_consistency(
    interchange, omm_system, omm_topology, pdb_file, print_detailed_info=False
):
    """
    Validate that the number of particles is consistent across different representations
    of the system.

    Args:
        interchange: OpenFF Interchange object
        omm_system: OpenMM System object
        omm_topology: OpenMM Topology object
    """
    # Get particle counts from different representations
    interchange_n_particles = len(interchange.positions)
    system_n_particles = omm_system.getNumParticles()
    topology_n_particles = omm_topology.getNumAtoms()

    # Hide just that exact message (case-sensitive)
    warnings.filterwarnings(
        "ignore",
        message=r"WARNING: duplicate atom",  # regex allowed
        module="openmm",  # safer than silencing everything
    )

    # Read the saved PDB file and get its particle count
    pdb = PDBFile(pdb_file)
    pdb_n_particles = pdb.topology.getNumAtoms()

    logger.info("System Consistency Check:")
    logger.info(f"Number of particles in Interchange: {interchange_n_particles}")
    logger.info(f"Number of particles in OpenMM System: {system_n_particles}")
    logger.info(f"Number of particles in OpenMM Topology: {topology_n_particles}")
    # print(f"Number of particles in saved PDB: {pdb_n_particles}")

    if print_detailed_info:
        # Additional diagnostic information
        logger.info("Detailed System Information:")
        # Check for NaN or undefined coordinates in interchange
        if hasattr(interchange, "positions"):
            nan_coords = np.isnan(interchange.positions).any(axis=1)
            if nan_coords.any():
                print(
                    f"Warning: Found {nan_coords.sum()} particles with NaN coordinates in interchange"
                )

        # Compare residue counts
        pdb_n_residues = pdb.topology.getNumResidues()
        omm_n_residues = omm_topology.getNumResidues()
        logger.info(f"Residue counts:")
        logger.info(f"OpenMM Topology: {omm_n_residues}")
        logger.info(f"PDB file: {pdb_n_residues}")

        # Compare chain counts
        pdb_n_chains = pdb.topology.getNumChains()
        omm_n_chains = omm_topology.getNumChains()
        logger.info(f"Chain counts:")
        logger.info(f"OpenMM Topology: {omm_n_chains}")
        logger.info(f"PDB file: {pdb_n_chains}")

        # Check for ligand (UNK residue)
        omm_unk_residues = [res for res in omm_topology.residues() if res.name == "UNK"]
        pdb_unk_residues = [res for res in pdb.topology.residues() if res.name == "UNK"]
        logger.info("Ligand (UNK) Information:")
        logger.info(
            f"Number of UNK residues in OpenMM Topology: {len(omm_unk_residues)}"
        )
        logger.info(f"Number of UNK residues in PDB: {len(pdb_unk_residues)}")

        if omm_unk_residues and pdb_unk_residues:
            # Compare atom counts in UNK residues
            omm_unk_atoms = sum(len(list(res.atoms())) for res in omm_unk_residues)
            pdb_unk_atoms = sum(len(list(res.atoms())) for res in pdb_unk_residues)
            logger.info(f"Number of atoms in UNK residues (OpenMM): {omm_unk_atoms}")
            logger.info(f"Number of atoms in UNK residues (PDB): {pdb_unk_atoms}")
        else:
            if not omm_unk_residues:
                logger.warning("Warning: No UNK residue found in OpenMM Topology!")
            if not pdb_unk_residues:
                logger.warning("Warning: No UNK residue found in PDB file!")

        # Count water molecules
        pdb_waters = sum(1 for res in pdb.topology.residues() if res.name == "HOH")
        omm_waters = sum(1 for res in omm_topology.residues() if res.name == "HOH")
        logger.info(f"Water molecule counts:")
        logger.info(f"OpenMM Topology: {omm_waters}")
        logger.info(f"PDB file: {pdb_waters}")

    # Check consistency across all representations
    if not (interchange_n_particles == system_n_particles == topology_n_particles):
        logger.error("Inconsistent number of particles found!")
        raise ValueError(
            "Inconsistent number of particles found!\n"
            f"Interchange: {interchange_n_particles}\n"
            f"OpenMM System: {system_n_particles}\n"
            f"OpenMM Topology: {topology_n_particles}\n"
            # f"Interchange PDB: {pdb_n_particles} {pdb_file}\n\n"
            "Please check the detailed system information above for potential causes."
        )
    else:
        logger.info("✅ Particle count is consistent across all representations")

    return True


# --------------------------------------------------------------------------
# Create OpenMM System and Topology Files
# --------------------------------------------------------------------------


@time_tracker
def save_openmm_system_topology(
    interchange,
    openff_interchange_path: str,
    openmm_topology_path: str,
    openmm_system_path: str,
    print_detailed_info,
):

    logger.info("Creating OpenMM System and Topology Files...")
    try:
        interchange.to_pdb(openff_interchange_path)

        omm_top = interchange.to_openmm_topology()
        with open(openmm_topology_path, "wb") as f:
            pickle.dump(omm_top, f)

        omm_system = interchange.to_openmm()

        with open(openmm_system_path, "w") as xml_file:
            xml_file.write(XmlSerializer.serialize(omm_system))

        # Perform system consistency validation
        validate_system_consistency(
            interchange,
            omm_system,
            omm_top,
            openff_interchange_path,
            print_detailed_info,
        )
        logger.info("✅ OpenMM System and Topology Files created successfully")
        logger.info(f"OpenFF Interchange: {openff_interchange_path}")
        logger.info(f"OpenMM Topology: {openmm_topology_path}")
        logger.info(f"OpenMM System: {openmm_system_path}")
    except Exception as e:
        logger.error(f"❌ Error creating OpenMM System and Topology Files: {e}")
        raise e

    return omm_system, omm_top


# --------------------------------------------------------------------------
# Parameterize OpenFF System
# --------------------------------------------------------------------------


@time_tracker
def parameterize_openff_system(openff_topology_path: str):
    """Create and parameterize OpenFF system from topology."""

    logger.info("Parameterizing OpenFF System...")
    try:
        with open(openff_topology_path, "r") as file:
            json_string = file.read()
        top = Topology.from_json(json_string)

        sage_ff14sb = ForceField(
            "openff-2.2.0.offxml",  # ← ligand parameters
            "ff14sb_off_impropers_0.0.3.offxml",
        )  # ← protein parameters

        interchange = sage_ff14sb.create_interchange(top)
        logger.info(
            f"✅ Forcefield parameterization complete. OpenFF interchange saved to: {openff_topology_path}"
        )
        return interchange
    except Exception as e:
        logger.error(f"❌ Error parameterizing OpenFF system: {e}")
        raise e


# --------------------------------------------------------------------------
# Create OpenFF Topology For Solvated System
# --------------------------------------------------------------------------


@time_tracker
def create_openff_topology(
    protein_topology: Topology,
    ligand_mol: Molecule,
    output_path: str,
    radius: Quantity = 2.5 * unit.angstrom,
    keep: Iterable[Molecule] = [],
) -> Topology:
    """Create complex topology by combining protein and ligand while handling clashes."""

    logger.info("Creating OpenFF Topology For Solvated System...")
    # If no ligand is provided, use the protein topology. No need to check for clashes with ligand.
    if ligand_mol is None:
        if check_for_large_molecules(protein_topology):
            with open(output_path, "w") as f:
                print(protein_topology.to_json(), file=f)
        logger.info(
            f"OpenFF Topology File saved to {output_path}. No ligand was provided, only protein was used."
        )
        return

    new_top_mols = []
    ligand_coordinates = ligand_mol.conformers[0][:, None, :]

    for molecule in protein_topology.molecules:
        # Keep molecules that are in the keep list
        if any(keep_mol.is_isomorphic_with(molecule) for keep_mol in keep):
            new_top_mols.append(molecule)
            continue

        if len(molecule.atoms) >= 50:
            new_top_mols.append(molecule)
            continue

        # For molecules not in the keep list, e.g. water molecules,
        # check if they are too close to the ligand. If they are, remove them.
        molecule_coordinates = molecule.conformers[0][None, :, :]
        diff_matrix = molecule_coordinates - ligand_coordinates
        working_unit = unit.nanometer
        distance_matrix = (
            np.linalg.norm(diff_matrix.m_as(working_unit), axis=-1) * working_unit
        )

        # If the molecule is not too close to the ligand, keep it.
        if distance_matrix.min() > radius:
            new_top_mols.append(molecule)
        else:
            logger.info(f"Removed {molecule.to_smiles()} molecule")

    # Add the ligand to the topology at the end
    new_top_mols.append(ligand_mol)

    # Create the new OpenFF Topology from the list of molecules
    new_top = Topology.from_molecules(new_top_mols)
    new_top.box_vectors = protein_topology.box_vectors

    # Check if there are any large molecules in the topology
    if check_for_large_molecules(new_top):
        with open(output_path, "w") as f:
            print(new_top.to_json(), file=f)
        logger.info(
            f"OpenFF Topology File saved to {output_path}. File includes protein and ligand."
        )
    else:
        raise ValueError(
            "❌ Protein was removed from topology. Topology saving failed."
        )


# --------------------------------------------------------------------------
# Main Function
# --------------------------------------------------------------------------


@time_tracker
def main(config, print_detailed_info=False):
    logger.info("========================================================")
    logger.info("🛠️ Forcefield Parameterization")
    logger.info("========================================================")

    # Check if ligand is provided in config.
    if config["path_ligand"] == "":
        ligand_mol = None
    else:
        ligand_mol = ligand_util.prepare_ligand_from_sdf(config["path_ligand"])

    protein_top_object = None  # Initialize first
    try:
        protein_top_object = Topology.from_pdb(config["path_protein_solvated"])
    except Exception as e:
        logger.error(f"❌ Error loading protein topology: {e}")
        return  # Exit the function if we can't load the protein topology

    create_openff_topology(
        protein_top_object, ligand_mol, config["path_openff_topology"]
    )

    interchange = parameterize_openff_system(config["path_openff_topology"])

    save_openmm_system_topology(
        interchange,
        config["path_openff_interchange"],
        config["path_openmm_topology"],
        config["path_openmm_system"],
        print_detailed_info,
    )


if __name__ == "__main__":
    main()
