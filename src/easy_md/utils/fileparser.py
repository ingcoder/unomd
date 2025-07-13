
# Standard library imports
import time
from functools import wraps
import subprocess


# Third-party imports
import MDAnalysis as mda
from rdkit import Chem
from openmm.app import *
from openmm import *
from openff.toolkit import ForceField, Topology
import openmm 
import parmed
import pickle
from openmm import XmlSerializer
from openmm.app import PDBFile, Simulation
from openff.toolkit import Topology
from openmm.unit import kelvin, picosecond, picoseconds
from openmm import(
    XmlSerializer,
    
)
from openff.toolkit import ForceField

# Custom imports & logging
from easy_md.utils import info_logger
import logging

logger = logging.getLogger(__name__)

def time_tracker(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"⏰ Function '{func.__name__}' took {execution_time:.4f} seconds to execute")
        return result
    return wrapper

def pqr_to_pdb(input_file_pqr, output_file_pdb):
    system = mda.Universe(input_file_pqr)
    system.atoms.write(output_file_pdb)

def pqr_to_pdb_obabel(input_file_pqr, output_file_pdb):
    command = ["obabel", input_file_pqr, "-O", output_file_pdb, "-xn"]
    try:
        subprocess.run(command, check=True)
        logger.info(f"Conversion successful! Output saved as {output_file_pdb}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during conversion: {e}")

def pdbqt_to_sdf(input_file, output_file, config):
    command = ["obabel", input_file, "-O", output_file, "-p", config['solv_pH'], "-h", "-m"]
    try:
        subprocess.run(command, check=True)
        logger.info(f"Conversion successful! Output saved as {output_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during conversion: {e}")

def pdbqt_to_pdb(input_file, output_file):
    command = ["obabel", input_file, "-O", output_file, "-m"]
    try:
        subprocess.run(command, check=True)
        logger.info(f"Conversion successful! Output saved as {output_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during conversion: {e}")
  
def prepare_ligand_for_md(input_file, output_file):
    """
    Ligand is needed to create the Protein-Ligand complex topology with OpenFF. 
    The OpenFF Molecule throws an error because it falsely detect free radicals in an aromatic ring. 
    1. Convert webina output.pdbqt to pdb with obabel. 
    2. This functions takes as input a pdb for a single position. 
    3. Postprocess it; fix valency and only add implicit hydrogens without changing the protonation used for docking. 
    4. Saves it as a SDF. 

    This approach ensures that the ligand SDF has all implict hydrogens, so it does not throw the free radical error when loaded with OpenFF Molecule. 
    It does not work converting webina.pdbqt directly to SDF with Openbabel. 
    """
    
    # Load PDB with existing hydrogens (from pH 6.2 preparation)
    mol = Chem.MolFromPDBFile(input_file, removeHs=False)

    # Sanitize the molecule (fix valency & bond perception)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL)

    # Add **only missing implicit hydrogens** (does not alter protonation state)
    mol = Chem.AddHs(mol, onlyOnAtoms=[a.GetIdx() for a in mol.GetAtoms() if a.GetTotalNumHs() == 0])

    # Save as an SDF
    writer = Chem.SDWriter(output_file)
    writer.write(mol)
    writer.close()
    print("✅ SDF saved with explicit and necessary implicit hydrogens!")

def interachange_top_to_amber_prmtop(config):
    top = Topology.from_json(open(config['path_openff_topology']).read())
    sage_ff14sb = ForceField(config['ff_small_molecule_openff'], config['ff_protein_openff'])
    interchange = sage_ff14sb.create_interchange(top)
    interchange.to_prmtop(config['path_amber_topology'], writer='internal')
    print(f"Done! File saved to {config['path_amber_topology']}")

def openmm_to_amber_topology(config):
    """
    Convert an OpenMM topology to an Amber topology using OpenFF Interchange.
    This ensures proper handling of force field parameters.
    
    Args:
        omm_system_xml_path (str): Path to the OpenMM system XML file
        omm_top_pkl_path (str): Path to the OpenMM topology pickle file
        off_topology_json_path (str): Path to the OpenFF topology JSON file
        output_prmtop_path (str): Path where the output Amber topology file will be saved
    """

    # Load the OpenMM system and topology
    off_top = load_openff_topology_from_json(config['path_openff_topology'])
    
    # Create an Interchange object from the OpenFF topology
    # Use the same force field combination you used to create the system
    sage_ff14sb = ForceField(config['ff_small_molecule_openff'], config['ff_protein_openff'])
    interchange = sage_ff14sb.create_interchange(off_top)
    
    # Convert to ParmEd Structure using the correct method
    parmed_structure = parmed.from_interchange(interchange)
    
    # Save as Amber prmtop file
    parmed_structure.save(config['path_amber_topology'], overwrite=True)
    print(f"✅ Successfully saved Amber topology to {config['path_amber_topology']}")
    
    return parmed_structure

def save_openmm_system_to_pdb(config):
    """ This function was used in "run_forcefield_parameterization.py:"  but is now deprecated and replaced by
    interchange.to_pdb(). The save_openmm_system_to_pdb function appears to be legacy code that's no longer used. It was likely replaced by the more integrated approach in save_openmm_system_topology which handles both system and topology conversion along with PDB file creation. I

    This function can still be used to load existing OpenMM system and save as PDB. 
    Save the initial positions from a openmm system file to PDB.
    
    Args:
        system_xml_path (str): Path to the OpenMM system XML file
        topology_pkl_path (str): Path to the OpenMM topology pickle file
        off_topology_json_path (str): Path to the OpenFF topology JSON file
        output_pdb_path (str): Path where the output PDB file will be saved
    """

    OMM_SYS_XML_INPUT = config['path_openmm_system']
    OMM_TOP_PKL_INPUT = config['path_openmm_topology']
    OFF_TOP_JSON_INPUT = config['path_openff_topology']
    PDB_OUTPUT = config['path_openmm_structure']
    
    omm_system = load_openmm_system_from_xml(OMM_SYS_XML_INPUT)
    omm_top = load_openmm_topology_from_pickle(OMM_TOP_PKL_INPUT)
    off_top = load_openff_topology_from_json(OFF_TOP_JSON_INPUT)        

    # 4) Integrator not used, but needed to create simulation object
    integrator = openmm.LangevinIntegrator(
        300 * kelvin, 
        1 / picosecond,
        0.002 * picoseconds,
    )

    # Combine the topology, system, integrator and initial positions into a simulation
    simulation = Simulation(omm_top, omm_system, integrator)
    simulation.context.setPositions(off_top.get_positions().to_openmm())

    # Write the positions to a PDB file
    with open(PDB_OUTPUT, 'w') as file:
        PDBFile.writeFile(simulation.topology, omm_top.get_positions().to_openmm(), file)

def load_openmm_system_from_xml(xml_file_path):
    with open(xml_file_path) as input:
        system = XmlSerializer.deserialize(input.read())
    return system

def load_openmm_topology_from_pickle(pickle_file_path):
    with open(pickle_file_path, 'rb') as f:
        return pickle.load(f)
    
def load_openff_topology_from_json(json_file_path):
    with open(json_file_path, 'r') as file:
        json_string = file.read()
    top = Topology.from_json(json_string) 
    return top


def load_files(system_xml_input, top_pkl_input, off_top_json_input=None):
    """Loads system and topology files required for the simulation."""
    with open(system_xml_input) as input:
        system = XmlSerializer.deserialize(input.read())
    
    with open(top_pkl_input, 'rb') as f:
        omm_top = pickle.load(f)

    if off_top_json_input:
        off_top = Topology.from_json(open(off_top_json_input).read())
    else:
        off_top = None

    return system, omm_top, off_top
