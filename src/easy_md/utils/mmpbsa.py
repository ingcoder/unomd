import parmed

def openmm_to_amber_topology(openmm_topology, amber_output_path):
    """
    Convert an OpenMM topology to an Amber topology.
    Assuming you have already created the Interchange object
    Generate OpenMM objects from the Interchange object
    """

    # Assuming you have already created the Interchange object
    # Generate OpenMM objects from the Interchange object
    # Load interchange 
    omm_system = interchange.to_openmm()
    omm_top = interchange.to_openmm_topology()

    # Convert to ParmEd Structure
    parmed_structure = parmed.openmm.load_topology(omm_top, omm_system)

    # Save as Amber prmtop file
    parmed_structure.save('output.prmtop', overwrite=True)
    pass

def check_residues_in_amber_topology(amber_topology_path):
    """
    Check if the residues in the Amber topology are correct.
    """
    # Print all unique residue names in the structure
    interchange.to_prmtop("out.prmtop") => ligand is missing in final prmtop
    interchange.to_inpcrd("out.inpcrd")

    amber_structure = parmed.load_file('out.prmtop', 'out.inpcrd')
    residue_names = set(residue.name for residue in amber_structure.residues)
    print("Residue names:", residue_names)


    # Check for water and ions
    water_residues = [res for res in amber_structure.residues if res.name in ['WAT', 'HOH']]
    sodium_ions = [res for res in amber_structure.residues if res.name == 'NA']
    chloride_ions = [res for res in amber_structure.residues if res.name == 'CL']

    print(f"Number of water molecules: {len(water_residues)}")
    print(f"Number of sodium ions: {len(sodium_ions)}")
    print(f"Number of chloride ions: {len(chloride_ions)}")
    pass

def openmm_to_amber_parameters(openmm_parameters):
    """
    Convert an OpenMM parameters to an Amber parameters.
    """
    pass