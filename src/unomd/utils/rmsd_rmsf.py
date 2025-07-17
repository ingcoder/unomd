import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.coordinates.DCD import DCDWriter 

# File Operations
def concat_traj(traj_files_input, top_file_input, concat_file_output):
    """Concatenate multiple trajectory files into a single DCD file."""
    # Check if all dcd files exist
    for file in traj_files_input:
        if not file.exists():
            raise FileNotFoundError(f"Directory {file} not found.")

    # Load universe from trajectory and topology
    u = mda.Universe(top_file_input, traj_files_input)

    try:         
        with DCDWriter(concat_file_output, n_atoms=u.atoms.n_atoms) as writer:
            for timestep in u.trajectory:
                writer.write(u.atoms)
        print(f"Combined trajectory saved to {concat_file_output}")
    except Exception as e:
        print(f"An error occurred: {e}")

def load_universe_and_reference(top_file_input, traj_file_input):
    """Load an MDAnalysis Universe and create a reference frame."""
    universe = mda.Universe(top_file_input, traj_file_input)
    # universe.trajectory[0]  # Ensure we are at the first frame
    # ref = universe.copy()  # Creates a copy of the universe at the first frame
    ref = mda.Universe(top_file_input)
    return universe, ref

# Atom Selection and Structure Analysis
def get_seg_ids(universe):
    """Get unique segment IDs for protein atoms in the universe."""
    protein_atoms = universe.select_atoms("protein")
    if len(protein_atoms) == 0:
        raise ValueError("No protein atoms found in the universe")
    return {atom.segid for atom in protein_atoms}

def align_universe_to_reference(universe, reference_frame, atoms_to_select, align_to_avg=False):
    """Align a trajectory to either a reference frame or average structure."""
    if align_to_avg:
        # Calculate and align to average structure (for RMSF)
        average = align.AverageStructure(universe, universe, 
                                       select=atoms_to_select, 
                                       ref_frame=0).run()
        ref = average.universe
        align.AlignTraj(universe, ref, 
                       select=atoms_to_select, 
                       in_memory=True).run()
    else:
        # Align to reference frame (for RMSD)
        align.AlignTraj(universe, reference_frame, 
                       select=atoms_to_select, 
                       in_memory=True).run()
    return universe

# RMSD Analysis
def calculate_rmsd(aligned_universe, atoms_to_select, ref):
    """Calculate RMSD between aligned trajectory and first frame."""
    # Create reference from first frame
    # aligned_universe.trajectory[0]
    # ref = aligned_universe.copy()

    # Calculate RMSD
    rmsd_calc = rms.RMSD(aligned_universe, ref, select=atoms_to_select)
    rmsd_calc.run()
    
    return pd.DataFrame(rmsd_calc.rmsd, 
                       columns=['Frame', 'Time (ps)', 'RMSD'])

def plot_rmsd_graph(rmsd_df):
    """Plot RMSD over time."""
    plt.figure(figsize=(10, 6))
    plt.plot(rmsd_df.index*0.2, rmsd_df['RMSD'], label='RMSD vs Time')
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD Over Time')
    plt.legend()
    plt.grid(True)
    plt.show()

def run_rmsd_analysis(top_file, traj_file, atoms_to_select):
    """
    Run complete RMSD analysis workflow.
    
    Parameters
    ----------
    top_file : str
        Path to topology file
    traj_file : str
        Path to trajectory file
    atoms_to_select : str
        MDAnalysis selection string for atoms to analyze
        
    Returns
    -------
    pandas.DataFrame
        DataFrame containing Frame, Time, and RMSD values
    """
    # Load trajectory and get reference frame
    universe, reference_frame = load_universe_and_reference(top_file, traj_file)
    
    # Align trajectory to reference frame
    aligned_universe = align_universe_to_reference(
        universe=universe,
        reference_frame=reference_frame,
        atoms_to_select=atoms_to_select,
        align_to_avg=False
    )
    
    # Calculate RMSD
    rmsd_df = calculate_rmsd(
        aligned_universe=aligned_universe,
        atoms_to_select=atoms_to_select,
        ref=reference_frame
    )
    
    return rmsd_df

# RMSF Analysis
def calculate_rmsf(aligned_universe, protein_segids):
    """Calculate RMSF for CA atoms in each protein segment."""
    rmsf_dict = {}
    
    for seg in protein_segids:
        # Select CA atoms for the current segment
        seg_atoms = aligned_universe.select_atoms(f'segid {seg} and protein and name CA')
        
        if len(seg_atoms) == 0:
            print(f"Warning: No CA atoms found for segment {seg}")
            continue
            
        # Calculate RMSF
        rmsf = rms.RMSF(seg_atoms).run()
        
        # Store results
        rmsf_dict[seg] = {
            'residues': seg_atoms.residues.resids,
            'rmsf': rmsf.results.rmsf
        }
    
    if not rmsf_dict:
        raise ValueError("No RMSF values could be calculated for any segment")
        
    return rmsf_dict

def plot_rmsf(rmsf_dict):
    """Plot RMSF for each protein segment."""
    fig, ax = plt.subplots(figsize=(10, 5))
    
    for key, value in rmsf_dict.items():
        resids = rmsf_dict[key]['residues']
        rmsf_vals = rmsf_dict[key]['rmsf']
        ax.plot(resids, rmsf_vals, label=f'Chain {key}')
    
    ax.set_xlabel('Unique Residue Number')
    ax.set_ylabel('RMSF')
    ax.set_title('RMSF for Each Chain')
    ax.legend()
    ax.grid(True)
    plt.show()

def run_rmsf_analysis(top_file, traj_file, atoms_to_select):
    """
    Run complete RMSF analysis workflow.
    
    Parameters
    ----------
    top_file : str
        Path to topology file
    traj_file : str
        Path to trajectory file
    atoms_to_select : str
        MDAnalysis selection string for atoms to align
        
    Returns
    -------
    dict
        Standardized results dictionary containing:
        - analysis_type: 'RMSF'
        - data: Dictionary of RMSF results per segment
        - selection: Atom selection string used
        - timestamp: When analysis was performed
        - metadata: Additional analysis information
    """
    # Load trajectory and get reference frame
    universe, reference_frame = load_universe_and_reference(top_file, traj_file)
    
    # Get protein segment IDs
    protein_segids = get_seg_ids(universe)
    
    # Align trajectory to average structure (required for RMSF)
    aligned_universe = align_universe_to_reference(
        universe=universe,
        reference_frame=reference_frame,
        atoms_to_select=atoms_to_select,
        align_to_avg=True # Always True for RMSF
    )
    
    # Calculate RMSF. Returns a dictionary of RMSF values for each protein segment
    return calculate_rmsf(
        aligned_universe=aligned_universe,
        protein_segids=protein_segids
    )