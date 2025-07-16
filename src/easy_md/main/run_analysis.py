#@title **Visualize RMSD/RMSF from MD trajectory**

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np
import warnings
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import os 

def run_analysis(config):

    warnings.filterwarnings('ignore')

    path_base = config.get("path_md_image")
    

    output_dir = config.get("analysis_dir")


    os.makedirs(output_dir, exist_ok=True)

    print("Loading trajectory files...")
    u = mda.Universe(f"{path_base}/output/emin.pdb", f"{path_base}/output/md_image_id_0.dcd")

    # === RMSD Calculation ===
    print("Calculating RMSD...")

    protein = u.select_atoms("protein")
    backbone = u.select_atoms("protein and backbone")

    rmsd_analysis = rms.RMSD(backbone, backbone, select="backbone", ref_frame=0)
    rmsd_analysis.run()

    ligand = u.select_atoms("resname EPE UNK LIG") # Adjust depending on ligand's name
    if len(ligand) > 0:
        rmsd_ligand = rms.RMSD(ligand, ligand, ref_frame=0)
        rmsd_ligand.run()
        print(f"Found ligand with {len(ligand)} atoms")

    print("Calculating RMSF...")

    align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()

    c_alphas = u.select_atoms("protein and name CA")
    positions = np.array([c_alphas.positions.copy() for ts in u.trajectory])
    avg_positions = positions.mean(axis=0)
    rmsf_values = np.sqrt(((positions - avg_positions)**2).sum(axis=2).mean(axis=0))
    residue_ids = c_alphas.resids

    # --- Plotly Visualization ---
    print("Generating interactive plots with Plotly...")

    fig = make_subplots(rows=2, cols=1,
                        subplot_titles=("Root Mean Square Deviation (RMSD)",
                                        "Root Mean Square Fluctuation (RMSF)"))

    # RMSD plot
    time = rmsd_analysis.rmsd[:, 1]
    rmsd_backbone = rmsd_analysis.rmsd[:, 2]

    fig.add_trace(go.Scatter(x=time, y=rmsd_backbone, mode='lines',
                            name='Protein Backbone',
                            line=dict(color='blue', width=2),
                            hovertemplate="Time: %{x:.2f} ps<br>RMSD: %{y:.2f} Å<extra></extra>"),
                row=1, col=1)

    if len(ligand) > 0:
        time_ligand = rmsd_ligand.rmsd[:, 1]
        rmsd_ligand_values = rmsd_ligand.rmsd[:, 2]
        fig.add_trace(go.Scatter(x=time_ligand, y=rmsd_ligand_values, mode='lines',
                                name='Ligand',
                                line=dict(color='red', width=2),
                                hovertemplate="Time: %{x:.2f} ps<br>RMSD: %{y:.2f} Å<extra></extra>"),
                    row=1, col=1)

    window = 10
    if len(time) > window:
        rmsd_smooth = np.convolve(rmsd_backbone, np.ones(window)/window, mode='valid')
        time_smooth = time[window//2:-window//2+1]
        fig.add_trace(go.Scatter(x=time_smooth, y=rmsd_smooth, mode='lines',
                                name=f'Protein Backbone (Smoothed, Window={window})',
                                line=dict(color='blue', width=3, dash='solid'),
                                opacity=0.5,
                                hovertemplate="Time: %{x:.2f} ps<br>RMSD (Smooth): %{y:.2f} Å<extra></extra>"),
                    row=1, col=1)

    mean_rmsd = np.mean(rmsd_backbone)
    std_rmsd = np.std(rmsd_backbone)

    fig.update_xaxes(title_text='Time (ps)', row=1, col=1)
    fig.update_yaxes(title_text='RMSD (Å)', row=1, col=1)
    fig.update_layout(hovermode="x unified")

    # RMSF plot
    fig.add_trace(go.Scatter(x=residue_ids, y=rmsf_values, mode='lines',
                            name='RMSF',
                            line=dict(color='green', width=2),
                            fill='tozeroy', fillcolor='rgba(0,255,0,0.2)',
                            hovertemplate="Residue: %{x}<br>RMSF: %{y:.2f} Å<extra></extra>"),
                row=2, col=1)

    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    threshold = mean_rmsf + 1.5 * std_rmsf

    fig.add_hline(y=threshold, line_dash="dash", line_color="red",
                annotation_text=f'Flexibility threshold ({threshold:.2f} Å)',
                annotation_position="top right", row=2, col=1, name='Threshold')

    high_flex_mask = rmsf_values > threshold
    if np.any(high_flex_mask):
        high_flex_res = residue_ids[high_flex_mask]
        high_flex_rmsf = rmsf_values[high_flex_mask]
        fig.add_trace(go.Scatter(x=high_flex_res, y=high_flex_rmsf, mode='markers',
                                name='Highly Flexible Residues',
                                marker=dict(color='red', size=8, symbol='circle'),
                                hovertemplate="Residue: %{x}<br>RMSF: %{y:.2f} Å<extra></extra>"),
                    row=2, col=1)

    fig.update_xaxes(title_text='Residue Number', row=2, col=1)
    fig.update_yaxes(title_text='RMSF (Å)', row=2, col=1)

    fig.update_layout(title_text='MD Trajectory Analysis: RMSD and RMSF', title_x=0.5,
                    height=800, showlegend=True)

    # Export as HTML
    html_plot_path = os.path.join(output_dir, "rmsd_rmsf_analysis_interactive.html")
    pio.write_html(fig, file=html_plot_path, auto_open=False)
    print(f"\n Interactive plot saved to: {html_plot_path}")

    # Display in Google Colab
    fig.show()

    print("\n" + "="*60)
    print("TRAJECTORY ANALYSIS SUMMARY")
    print("="*60)
    print(f"Simulation length: {time[-1]:.2f} ps")
    print(f"Number of frames: {len(u.trajectory)}")
    print(f"Time step between frames: {time[1]-time[0]:.4f} ps")

    print(f"\n PROTEIN BACKBONE RMSD:")
    print(f"   Average: {mean_rmsd:.3f} ± {std_rmsd:.3f} Å")
    print(f"   Range: {np.min(rmsd_backbone):.3f} - {np.max(rmsd_backbone):.3f} Å")
    print(f"   Final RMSD: {rmsd_backbone[-1]:.3f} Å")

    if len(ligand) > 0:
        print(f"\n LIGAND RMSD:")
        print(f"   Average: {np.mean(rmsd_ligand_values):.3f} ± {np.std(rmsd_ligand_values):.3f} Å")
        print(f"   Range: {np.min(rmsd_ligand_values):.3f} - {np.max(rmsd_ligand_values):.3f} Å")

    print(f"\n RESIDUE FLEXIBILITY (RMSF):")
    print(f"   Average: {mean_rmsf:.3f} ± {std_rmsf:.3f} Å")
    print(f"   Most flexible: Residue {residue_ids[np.argmax(rmsf_values)]} ({np.max(rmsf_values):.3f} Å)")
    print(f"   Most stable: Residue {residue_ids[np.argmin(rmsf_values)]} ({np.min(rmsf_values):.3f} Å)")

    if np.any(high_flex_mask):
        print(f"\n   Highly flexible residues (>{threshold:.2f} Å):")
        for res, flex in zip(high_flex_res, high_flex_rmsf):
            print(f"     - Residue {res}: {flex:.3f} Å")

    # Save data files to the new output_dir
    np.savetxt(os.path.join(output_dir, "rmsd_data.txt"),
            np.column_stack([time, rmsd_backbone]),
            header="Time(ps) RMSD(Angstrom)",
            fmt='%.4f')

    np.savetxt(os.path.join(output_dir, "rmsf_data.txt"),
            np.column_stack([residue_ids, rmsf_values]),
            header="Residue RMSF(Angstrom)",
            fmt='%d %.4f')

    print(f"\nData saved as: {os.path.join(output_dir, 'rmsd_data.txt')} and {os.path.join(output_dir, 'rmsf_data.txt')}")