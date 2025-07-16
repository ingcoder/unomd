#@title **Visualize MD trajectory**

!pip install py3Dmol -q
!pip install MDAnalysis -q
import py3Dmol
import MDAnalysis as mda
import warnings
import os

def run_analysis(config):

    warnings.filterwarnings('ignore', message="Found no information for attr: 'formalcharges'")
    warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')
    warnings.filterwarnings('ignore', category=DeprecationWarning, module='MDAnalysis.coordinates.DCD')

    path_base = "/content/easy-md/example"
    output_dir = f"{path_base}/output/analysis"

    emin_structure = f"{output_dir}/emin.pdb"
    md_trajectory = f"{output_dir}/md_image_id_0.dcd"

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

        html_output_path = f"{output_dir}/trajectory_animation_3D.html"
        html = view._make_html()
        with open(html_output_path, 'w') as f:
            f.write(html)

        print(f"Location: {html_output_path}")
        print(f"The HTML file contains the full animated trajectory")
        print(f"Total frames in animation: {len(u.trajectory)}")