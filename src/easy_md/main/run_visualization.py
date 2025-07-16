# Visualize MD trajectory
# Required dependencies: py3Dmol, MDAnalysis
# Install with: pip install py3Dmol MDAnalysis
import py3Dmol
import MDAnalysis as mda
import warnings
import os
import glob
import logging

logger = logging.getLogger(__name__)

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

def run_analysis(config):

    warnings.filterwarnings('ignore', message="Found no information for attr: 'formalcharges'")
    warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')
    warnings.filterwarnings('ignore', category=DeprecationWarning, module='MDAnalysis.coordinates.DCD')

    output_dir = config.get("analysis_dir")

    emin_structure = config.get("path_emin_structure")
    md_trajectory = get_last_final_dcd(config)

    if not os.path.exists(emin_structure) or not os.path.exists(md_trajectory):
        print("❌ MD simulation files not found. Please run the quickrun() command first.")
    else:
        print("Loading trajectory files...")

        u = mda.Universe(emin_structure, md_trajectory)
        system = u.select_atoms('not (resname HOH CL NA SOL WAT TIP3P)')

        print(f"Creating animation with {len(u.trajectory)} frames...")

        temp_traj = '/tmp/trajectory_for_viz.pdb'
        system.write(temp_traj, frames='all')

        view = py3Dmol.view(width=800, height=600)
        with open(temp_traj) as f:
            view.addModelsAsFrames(f.read())

        view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})
        view.setStyle({'model': -1, 'resn': ['UNK', 'EPE']}, {'stick': {'radius': 0.3, 'color': 'red'}})
        view.animate({'loop': 'forward', 'interval': 75})
        view.zoomTo()
        print("\n=== 3D Trajectory Animation ===")
        view.show()

        html_output_path = os.path.join(output_dir, "trajectory_animation_3D.html")
        html = view._make_html()
        with open(html_output_path, 'w') as f:
            f.write(html)

        print(f"Location: {html_output_path}")
        print(f"The HTML file contains the full animated trajectory")
        print(f"Total frames in animation: {len(u.trajectory)}")