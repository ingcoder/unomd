"""
This script provides utilities for analyzing and validating force field parameters in molecular systems.
It focuses on examining nonbonded interactions, detecting atomic clashes between protein and ligand 
molecules, and checking simulation system configurations.

In simpler terms, this script helps:
- Check if atoms in a molecular model are positioned correctly without unrealistic overlaps
- Extract and analyze parameters that describe how atoms interact with each other
- Validate settings for computer simulations of molecules

These validations are important before running computationally expensive molecular simulations
to ensure the results will be reliable and scientifically sound.
"""

from openmm import NonbondedForce, XmlSerializer
import numpy as np
from Bio.PDB import PDBParser
from openmm import unit
from openmm import unit
import warnings
warnings.filterwarnings("ignore")


def get_sigma_epsilon(nonbonded_force, particle_index):
    """
    Extract and correct Lennard-Jones parameters for a particle.
    
    Purpose:
        Gets parameters that describe how atoms interact when they're not bonded to each other.
    
    What are Lennard-Jones parameters?
        - Sigma: Represents the size of an atom (distance where interaction energy is zero)
        - Epsilon: Represents how strongly atoms attract each other
        These parameters model how atoms repel at close distances but attract at moderate distances.
    
    Why adjustment is needed:
        OpenMM uses a slightly different definition of sigma than the standard Lennard-Jones
        equation. We adjust sigma by multiplying by 2^(1/6) to convert between these 
        conventions. Without this correction, calculations of atom overlaps would be inaccurate.
        
    Args:
        nonbonded_force: OpenMM NonbondedForce object
        particle_index: Index of the particle
        
    Returns:
        sigma_corrected: Corrected sigma value (atom size parameter)
        epsilon: Epsilon value (attraction strength parameter)
    """
    charge, sigma, epsilon = nonbonded_force.getParticleParameters(particle_index)
    sigma_corrected = sigma * (2 ** (1/6)) # Adjust sigma to match the standard Lennard-Jones sigma
    return sigma_corrected, epsilon


def find_clashes(protein_atoms, ligand_coords, nonbonded_force):
    """
    Identify atomic clashes between protein and ligand atoms using Lennard-Jones parameters.
    The sigma value, as defined in the standard Lennard-Jones formula, represents the distance at which the interaction 
    energy between two non-bonded atoms is zero - essentially defining the "effective size" of an atom. 
    This corresponds to the van der Waals radius, which describes how close atoms can approach each other 
    before strong repulsion occurs.

    
    Purpose:
        Finds places where atoms are too close to each other, which isn't physically realistic.
    
    Why useful:
        In real molecules, atoms can't overlap significantly due to electron repulsion.
        Finding these "clashes" helps identify problems in the molecular model that need
        to be fixed before simulation. Think of it like finding places where two billiard
        balls are trying to occupy the same space - that can't happen in reality.
    
    Args:
        protein_atoms: List of protein atoms
        ligand_coords: Coordinates of ligand atoms
        nonbonded_force: OpenMM NonbondedForce object
        
    Returns:
        clash_pairs: List of details about detected clashes
    """
    # Ensure coordinates are in OpenMM units (nanometers typically)
    protein_coords = np.array([atom.coord for atom in protein_atoms]) * unit.angstroms  # if coords are in angstroms
    clash_pairs = []

    # Iterate over all ligand atoms and check if clashes with protein atoms, based on sigma
    for lig_idx, lig_coord in enumerate(ligand_coords):
        lig_coord = np.array(lig_coord) * unit.angstroms  # Convert ligand coords to proper OpenMM units
        lig_sigma, lig_epsilon = get_sigma_epsilon(nonbonded_force, lig_idx + 188832)  # Adjust index appropriately
        # lig_sigma *= unit.nanometers  # Ensure sigma is in nanometers if needed
    
        for prot_idx, prot_coord in enumerate(protein_coords):
            prot_sigma, prot_epsilon = get_sigma_epsilon(nonbonded_force, prot_idx)
            # prot_sigma *= unit.nanometers  # Ensure sigma is in nanometers if needed
            if lig_idx == 0:
                if prot_idx == 0:
                    radius_sum = lig_sigma + prot_sigma
                    distance = np.linalg.norm(prot_coord - lig_coord) * unit.nanometers
                    distance /=10
            radius_sum = (lig_sigma + prot_sigma) / 2
            distance = np.linalg.norm(prot_coord - lig_coord) * unit.nanometers # No unit conversion needed if both are already in nanometers
            distance /=10
            
            if distance <= radius_sum:
                overlap = (radius_sum - distance)
                if overlap > 0.1 * unit.nanometers:
                    print(f"lig_sigma {lig_sigma}")
                    print(f"prot_sigma {prot_sigma}")
                    clash_pairs.append(
                        (f"Protein Atom Index: {prot_idx}",
                        f"Ligand Atom Index: {lig_idx}",
                        f"Distance: {distance * unit.nanometers} Å",
                        f"Radius Sum: {radius_sum  * unit.nanometers} Å",
                        f"Overlap: {(radius_sum - distance) * unit.nanometers} Å")
                    )
        print(f"Number of clashes: {len(clash_pairs)}")
        return clash_pairs

def extract_protein_and_ligand_atoms(pdb_file, protein_chain_ids, ligand_resname):
    """
    Extract protein and ligand atoms from a PDB file.
    
    Purpose:
        Separates the protein parts from the drug/ligand parts in a molecular structure file.
    
    Why useful:
        It's like sorting the pieces of a complex puzzle into groups. This separation
        allows us to analyze how the drug interacts with the protein, making it easier
        to focus on the important interactions without getting lost in all the details
        of the entire molecular system.
    
    Args:
        pdb_file: Path to the PDB file (a standard file format for 3D molecular structures)
        protein_chain_ids: List of chain IDs to extract protein atoms from
        ligand_resname: Name that identifies the ligand in the file
        
    Returns:
        protein_atoms: List of protein atoms
        ligand_coords: Coordinates of ligand atoms
        ligand_elements: Chemical elements of ligand atoms
    """
    parser = PDBParser()
    structure = parser.get_structure('ID', pdb_file)
    # Select the first model
    model = structure[0]
    # Extracting protein atoms from specified chains
    protein_atoms = [atom for atom in model.get_atoms() 
                     if atom.get_parent().get_parent().id in protein_chain_ids]
    # Extracting ligand atoms (assuming the ligand is identified by a specific residue name)
    ligand_atoms = [atom for atom in model.get_atoms() 
                    if atom.get_parent().resname == ligand_resname]
    # Extract ligand coordinates and elements
    ligand_coords = [atom.get_coord() for atom in ligand_atoms]
    ligand_elements = [atom.element for atom in ligand_atoms]
    print(f"Protein length: {len(protein_atoms)}")

    return protein_atoms, ligand_coords, ligand_elements


# Check if the system uses periodic boundary conditions
def check_periodic_boundaries(system):
    """
    Check and print periodic boundary conditions of an OpenMM system.
    
    Purpose:
        Examines how the simulation handles the edges of the simulation box.
    
    Why useful:
        In simulations, molecules are placed in a box. Periodic boundaries act like
        portals - when a molecule leaves one side of the box, it reappears on the
        opposite side. This prevents unwanted edge effects and mimics a larger system.
        Checking these settings ensures the simulation won't have artificial boundary
        problems, especially for simulations in water or other solvents.
    
    Args:
        system: OpenMM System object
    """
    box_vectors = system.getDefaultPeriodicBoxVectors()
    if box_vectors is None:
        print("No periodic boundary conditions are set.")
    else:
        print("Periodic boundary conditions are set.")
        print("Box Vectors:")
        for vec in box_vectors:
            print(vec)

def run_molecule_analysis():
    """
    The system's force objects contain the mathematical rules and parameters that govern molecular interactions in the 
    simulation. 
    
    The NonbondedForce object specifically contains:
    - Electrostatic parameters:
        Atomic charges for each particle
        Method for handling long-range electrostatic interactions (PME - Particle Mesh Ewald)
    - Van der Waals parameters:
        Sigma values (atomic size parameters)
        Epsilon values (interaction strength parameters)
    - Interaction cutoffs:
        Distance cutoffs for calculations (beyond which interactions are ignored)
        Method for handling interactions (shown in "Current Nonbonded Method")
    - Particle information:
        Total count of particles in the system
        Individual parameter sets for each atom
    """
    
    # Set Variables
    protein_chain_ids = ['A', 'B']  # Chains from which to extract protein atoms
    ligand_resname = 'UNK'          # Residue name for the ligand
    openmm_system_filepath = '../files/ATP5IF1/complex/complex_if1_dimer_ua_ph6.2_faspr_solvated_system3.xml'
    pdb_filepath = '../files/ATP5IF1/complex/complex_if1_dimer_ua_ph6.2_faspr_solvated_emin.pdb' 

    with open(openmm_system_filepath) as input:
        system = XmlSerializer.deserialize(input.read())

    # Access NonbondedForce parametersfrom the system
    nonbonded_force = None
    for force in system.getForces():
        force_type = type(force).__name__
        print(f"\n{force_type} Parameters:")
        if isinstance(force, NonbondedForce):
            nonbonded_force = force
            method = force.getNonbondedMethod()
            print(force.PME)
            print("Current Nonbonded Method:", type(method))
            print("Cutoff distance for nonbonded interactions (nm):", force.getCutoffDistance())
            break
    if nonbonded_force is None:
        raise ValueError("NonbondedForce not found in the system")
    
    print(f"Total number of particles in the nonbonded_force: {nonbonded_force.getNumParticles()}")
    protein_atoms, ligand_coords, _ = extract_protein_and_ligand_atoms(pdb_filepath, protein_chain_ids, ligand_resname)
    find_clashes(protein_atoms, ligand_coords, nonbonded_force)
    check_periodic_boundaries(system)
    

if __name__=="__main__":
    run_molecule_analysis()
