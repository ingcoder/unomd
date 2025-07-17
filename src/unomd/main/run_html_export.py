"""
Comprehensive HTML export and analysis suite for MD trajectories.

This module provides a complete toolkit for analyzing and visualizing molecular dynamics 
trajectories with interactive HTML reports. It includes:

- RMSD (Root Mean Square Deviation) analysis for protein backbone and ligands
- RMSF (Root Mean Square Fluctuation) analysis for residue flexibility
- Interactive Plotly visualizations with statistical summaries
- 3D trajectory animation using py3Dmol with customizable styling
- Comprehensive HTML analysis reports with summary tables
- Data export functionality for further analysis
- Automated identification of highly flexible residues
- Statistical analysis including convergence assessment

The module generates publication-ready HTML reports that can be easily shared
and viewed in any web browser, making it ideal for collaborative research
and presentation purposes.
"""

# Standard library imports
import os
import warnings
from pathlib import Path
import yaml
import glob

# Third-party imports
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import plotly.express as px
import py3Dmol

# Custom imports
from unomd.utils.fileparser import time_tracker
from unomd.utils import info_logger
import logging


logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------


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

def _calculate_rmsd(u, ligand_selection="resname EPE UNK LIG"):
    """Calculate RMSD for protein backbone and ligand."""
    # Defensive check to ensure universe is properly loaded
    if not hasattr(u, 'atoms') or len(u.atoms) == 0:
        raise ValueError("Universe object is not properly loaded or has no atoms")
    
    backbone = u.select_atoms("protein and backbone")
    if len(backbone) == 0:
        raise ValueError("No protein backbone atoms found in the structure")
    
    rmsd_analysis = rms.RMSD(backbone, backbone, select="backbone", ref_frame=0)
    rmsd_analysis.run()
    
    ligand = u.select_atoms(ligand_selection)
    rmsd_ligand = None
    if len(ligand) > 0:
        rmsd_ligand = rms.RMSD(ligand, ligand, ref_frame=0)
        rmsd_ligand.run()
        logger.info(f"Found ligand with {len(ligand)} atoms")
    
    return rmsd_analysis, rmsd_ligand


def _calculate_rmsf(u):
    """Calculate RMSF for protein C-alpha atoms."""
    # Defensive check to ensure universe is properly loaded
    if not hasattr(u, 'atoms') or len(u.atoms) == 0:
        raise ValueError("Universe object is not properly loaded or has no atoms")
    
    # Check if there are C-alpha atoms available
    c_alphas = u.select_atoms("protein and name CA")
    if len(c_alphas) == 0:
        raise ValueError("No protein C-alpha atoms found in the structure")
    
    align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
    
    positions = np.array([c_alphas.positions.copy() for ts in u.trajectory])
    avg_positions = positions.mean(axis=0)
    rmsf_values = np.sqrt(((positions - avg_positions)**2).sum(axis=2).mean(axis=0))
    
    return rmsf_values, c_alphas.resids


def _create_interactive_plot(rmsd_analysis, rmsd_ligand, rmsf_values, residue_ids):
    """Create enhanced interactive Plotly visualization."""
    # Create subplots with improved spacing
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "Root Mean Square Deviation (RMSD) vs Time", 
            "RMSD Statistics", 
            "Root Mean Square Fluctuation (RMSF) vs Residue",
            "RMSF Distribution"
        ),
        specs=[[{"colspan": 2}, None], [{"colspan": 2}, None]],
        vertical_spacing=0.12,
        horizontal_spacing=0.05
    )
    
    # Extract data
    time = rmsd_analysis.rmsd[:, 1]
    rmsd_backbone = rmsd_analysis.rmsd[:, 2]
    
    # Calculate statistics
    mean_rmsd = np.mean(rmsd_backbone)
    std_rmsd = np.std(rmsd_backbone)
    
    # Enhanced RMSD plot
    fig.add_trace(
        go.Scatter(
            x=time, y=rmsd_backbone, 
            mode='lines+markers', 
            name='Protein Backbone',
            line=dict(color='#2E86AB', width=3),
            marker=dict(size=4, color='#2E86AB', opacity=0.7),
            hovertemplate="<b>Time:</b> %{x:.3f} ps<br><b>RMSD:</b> %{y:.3f} Å<extra></extra>"
        ),
        row=1, col=1
    )
    
    # Add mean line
    fig.add_hline(
        y=mean_rmsd, 
        line_dash="dash", 
        line_color="rgba(46, 134, 171, 0.8)",
        annotation_text=f'Mean: {mean_rmsd:.3f} Å',
        annotation_position="top left", 
        row=1, col=1
    )
    
    # Add ligand RMSD if available
    if rmsd_ligand is not None:
        time_ligand = rmsd_ligand.rmsd[:, 1]
        rmsd_ligand_values = rmsd_ligand.rmsd[:, 2]
        mean_ligand = np.mean(rmsd_ligand_values)
        
        fig.add_trace(
            go.Scatter(
                x=time_ligand, y=rmsd_ligand_values, 
                mode='lines+markers', 
                name='Ligand',
                line=dict(color='#A23B72', width=3),
                marker=dict(size=4, color='#A23B72', opacity=0.7),
                hovertemplate="<b>Time:</b> %{x:.3f} ps<br><b>RMSD:</b> %{y:.3f} Å<extra></extra>"
            ),
            row=1, col=1
        )
        
        fig.add_hline(
            y=mean_ligand, 
            line_dash="dash", 
            line_color="rgba(162, 59, 114, 0.8)",
            annotation_text=f'Ligand Mean: {mean_ligand:.3f} Å',
            annotation_position="top right", 
            row=1, col=1
        )
    
    # Enhanced RMSF plot
    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    threshold = mean_rmsf + 1.5 * std_rmsf
    
    fig.add_trace(
        go.Scatter(
            x=residue_ids, y=rmsf_values, 
            mode='lines+markers', 
            name='RMSF',
            line=dict(color='#F18F01', width=2), 
            marker=dict(size=3, color='#F18F01', opacity=0.8),
            fill='tozeroy', 
            fillcolor='rgba(241, 143, 1, 0.3)',
            hovertemplate="<b>Residue:</b> %{x}<br><b>RMSF:</b> %{y:.3f} Å<extra></extra>"
        ),
        row=2, col=1
    )
    
    # Add flexibility threshold
    fig.add_hline(
        y=threshold, 
        line_dash="dash", 
        line_color="#C73E1D",
        annotation_text=f'Flexibility Threshold: {threshold:.3f} Å',
        annotation_position="top right", 
        row=2, col=1
    )
    
    # Add mean line for RMSF
    fig.add_hline(
        y=mean_rmsf, 
        line_dash="dot", 
        line_color="rgba(241, 143, 1, 0.8)",
        annotation_text=f'Mean: {mean_rmsf:.3f} Å',
        annotation_position="bottom right", 
        row=2, col=1
    )
    
    # Highlight highly flexible residues
    high_flex_mask = rmsf_values > threshold
    if np.any(high_flex_mask):
        high_flex_res = residue_ids[high_flex_mask]
        high_flex_rmsf = rmsf_values[high_flex_mask]
        fig.add_trace(
            go.Scatter(
                x=high_flex_res, y=high_flex_rmsf, 
                mode='markers',
                name='Highly Flexible',
                marker=dict(color='#C73E1D', size=8, symbol='diamond'),
                hovertemplate="<b>Residue:</b> %{x}<br><b>RMSF:</b> %{y:.3f} Å<extra></extra>"
            ),
            row=2, col=1
        )
    
    # Update layout with enhanced styling
    fig.update_xaxes(title_text='Time (ps)', row=1, col=1, showgrid=True, gridwidth=1, gridcolor='rgba(128,128,128,0.2)')
    fig.update_yaxes(title_text='RMSD (Å)', row=1, col=1, showgrid=True, gridwidth=1, gridcolor='rgba(128,128,128,0.2)')
    fig.update_xaxes(title_text='Residue Number', row=2, col=1, showgrid=True, gridwidth=1, gridcolor='rgba(128,128,128,0.2)')
    fig.update_yaxes(title_text='RMSF (Å)', row=2, col=1, showgrid=True, gridwidth=1, gridcolor='rgba(128,128,128,0.2)')
    
    fig.update_layout(
        title={
            'text': 'MD Trajectory Analysis: RMSD and RMSF',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 20, 'color': '#2E86AB'}
        },
        height=900, 
        showlegend=True, 
        hovermode="x unified",
        plot_bgcolor='rgba(248,248,248,0.8)',
        paper_bgcolor='white',
        font=dict(size=12),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    return fig


def _create_analysis_summary_table(u, rmsd_analysis, rmsd_ligand, rmsf_values, residue_ids):
    """Create a comprehensive analysis summary table."""
    time = rmsd_analysis.rmsd[:, 1]
    rmsd_backbone = rmsd_analysis.rmsd[:, 2]
    
    # Calculate statistics
    mean_rmsd = np.mean(rmsd_backbone)
    std_rmsd = np.std(rmsd_backbone)
    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    threshold = mean_rmsf + 1.5 * std_rmsf
    
    # Create summary data
    summary_data = {
        'Simulation Parameters': [
            ['Simulation Length', f'{time[-1]:.2f} ps'],
            ['Number of Frames', f'{len(u.trajectory)}'],
            ['Time Step', f'{time[1]-time[0]:.4f} ps'],
            ['Total Atoms', f'{len(u.atoms)}'],
            ['Protein Atoms', f'{len(u.select_atoms("protein"))}']
        ],
        'RMSD Statistics (Å)': [
            ['Backbone Mean ± STD', f'{mean_rmsd:.3f} ± {std_rmsd:.3f}'],
            ['Backbone Range', f'{np.min(rmsd_backbone):.3f} - {np.max(rmsd_backbone):.3f}'],
            ['Final RMSD', f'{rmsd_backbone[-1]:.3f}'],
            ['Convergence Status', 'Converged' if rmsd_backbone[-1] < mean_rmsd + std_rmsd else 'Still Equilibrating']
        ],
        'RMSF Statistics (Å)': [
            ['Mean ± STD', f'{mean_rmsf:.3f} ± {std_rmsf:.3f}'],
            ['Range', f'{np.min(rmsf_values):.3f} - {np.max(rmsf_values):.3f}'],
            ['Most Flexible', f'Residue {residue_ids[np.argmax(rmsf_values)]} ({np.max(rmsf_values):.3f})'],
            ['Most Stable', f'Residue {residue_ids[np.argmin(rmsf_values)]} ({np.min(rmsf_values):.3f})'],
            ['Flexibility Threshold', f'{threshold:.3f}']
        ]
    }
    
    # Add ligand statistics if available
    if rmsd_ligand is not None:
        rmsd_ligand_values = rmsd_ligand.rmsd[:, 2]
        summary_data['Ligand RMSD (Å)'] = [
            ['Mean ± STD', f'{np.mean(rmsd_ligand_values):.3f} ± {np.std(rmsd_ligand_values):.3f}'],
            ['Range', f'{np.min(rmsd_ligand_values):.3f} - {np.max(rmsd_ligand_values):.3f}'],
            ['Final RMSD', f'{rmsd_ligand_values[-1]:.3f}'],
            ['Stability', 'Stable' if np.std(rmsd_ligand_values) < 0.5 else 'Mobile']
        ]
    
    # Create HTML table
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>MD Analysis Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
            .container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
            h1 { color: #2E86AB; text-align: center; margin-bottom: 30px; }
            .table-container { margin-bottom: 30px; }
            .section-title { color: #2E86AB; font-size: 18px; font-weight: bold; margin-bottom: 10px; border-bottom: 2px solid #2E86AB; padding-bottom: 5px; }
            table { width: 100%; border-collapse: collapse; margin-bottom: 20px; }
            th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
            th { background-color: #2E86AB; color: white; font-weight: bold; }
            tr:nth-child(even) { background-color: #f9f9f9; }
            tr:hover { background-color: #f5f5f5; }
            .metric { font-weight: bold; color: #333; }
            .value { color: #666; }
            .flexible-residues { background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 15px; border-radius: 5px; margin-top: 10px; }
            .flexible-residues h3 { color: #C73E1D; margin-top: 0; }
            .residue-list { display: flex; flex-wrap: wrap; gap: 10px; }
            .residue-item { background-color: #C73E1D; color: white; padding: 5px 10px; border-radius: 15px; font-size: 12px; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>MD Trajectory Analysis Summary</h1>
    """
    
    # Add tables for each section
    for section_title, data in summary_data.items():
        html_content += f"""
            <div class="table-container">
                <div class="section-title">{section_title}</div>
                <table>
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th>Value</th>
                        </tr>
                    </thead>
                    <tbody>
        """
        
        for metric, value in data:
            html_content += f"""
                        <tr>
                            <td class="metric">{metric}</td>
                            <td class="value">{value}</td>
                        </tr>
            """
        
        html_content += """
                    </tbody>
                </table>
            </div>
        """
    
    # Add highly flexible residues section
    high_flex_mask = rmsf_values > threshold
    if np.any(high_flex_mask):
        high_flex_res = residue_ids[high_flex_mask]
        high_flex_rmsf = rmsf_values[high_flex_mask]
        
        html_content += f"""
            <div class="flexible-residues">
                <h3>Highly Flexible Residues (>{threshold:.3f} Å)</h3>
                <div class="residue-list">
        """
        
        for res, flex in zip(high_flex_res, high_flex_rmsf):
            html_content += f'<div class="residue-item">Res {res}: {flex:.3f} Å</div>'
        
        html_content += """
                </div>
            </div>
        """
    
    html_content += """
        </div>
    </body>
    </html>
    """
    
    return html_content


def _save_analysis_data(config, rmsd_analysis, rmsf_values, residue_ids, rmsd_ligand=None):
    """Save analysis data to text files."""

    
    # Save RMSD data
    time = rmsd_analysis.rmsd[:, 1]
    rmsd_backbone = rmsd_analysis.rmsd[:, 2]
    np.savetxt(
        config.get("path_rmsd_output_text"),
        np.column_stack([time, rmsd_backbone]),
        header="Time(ps) RMSD(Angstrom)",
        fmt='%.4f'
    )
    
    # Save RMSF data
    np.savetxt(
        config.get("path_rmsf_output_text"),
        np.column_stack([residue_ids, rmsf_values]),
        header="Residue RMSF(Angstrom)",
        fmt='%d %.4f'
    )
    
    # Save ligand RMSD if available
    if rmsd_ligand is not None:
        time_ligand = rmsd_ligand.rmsd[:, 1]
        rmsd_ligand_values = rmsd_ligand.rmsd[:, 2]
        np.savetxt(
            config.get("path_rmsd_ligand_output_text"),
            np.column_stack([time_ligand, rmsd_ligand_values]),
            header="Time(ps) RMSD_Ligand(Angstrom)",
            fmt='%.4f'
        )


def _log_analysis_summary(u, rmsd_analysis, rmsd_ligand, rmsf_values, residue_ids):
    """Log comprehensive analysis summary."""
    time = rmsd_analysis.rmsd[:, 1]
    rmsd_backbone = rmsd_analysis.rmsd[:, 2]
    
    logger.info("=" * 60)
    logger.info("TRAJECTORY ANALYSIS SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Simulation length: {time[-1]:.2f} ps")
    logger.info(f"Number of frames: {len(u.trajectory)}")
    logger.info(f"Time step between frames: {time[1]-time[0]:.4f} ps")
    
    mean_rmsd = np.mean(rmsd_backbone)
    std_rmsd = np.std(rmsd_backbone)
    logger.info(f"PROTEIN BACKBONE RMSD:")
    logger.info(f"   Average: {mean_rmsd:.3f} ± {std_rmsd:.3f} Å")
    logger.info(f"   Range: {np.min(rmsd_backbone):.3f} - {np.max(rmsd_backbone):.3f} Å")
    logger.info(f"   Final RMSD: {rmsd_backbone[-1]:.3f} Å")
    
    if rmsd_ligand is not None:
        rmsd_ligand_values = rmsd_ligand.rmsd[:, 2]
        logger.info(f"LIGAND RMSD:")
        logger.info(f"   Average: {np.mean(rmsd_ligand_values):.3f} ± {np.std(rmsd_ligand_values):.3f} Å")
        logger.info(f"   Range: {np.min(rmsd_ligand_values):.3f} - {np.max(rmsd_ligand_values):.3f} Å")
    
    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    threshold = mean_rmsf + 1.5 * std_rmsf
    
    logger.info(f"RESIDUE FLEXIBILITY (RMSF):")
    logger.info(f"   Average: {mean_rmsf:.3f} ± {std_rmsf:.3f} Å")
    logger.info(f"   Most flexible: Residue {residue_ids[np.argmax(rmsf_values)]} ({np.max(rmsf_values):.3f} Å)")
    logger.info(f"   Most stable: Residue {residue_ids[np.argmin(rmsf_values)]} ({np.min(rmsf_values):.3f} Å)")
    
    high_flex_mask = rmsf_values > threshold
    if np.any(high_flex_mask):
        high_flex_res = residue_ids[high_flex_mask]
        high_flex_rmsf = rmsf_values[high_flex_mask]
        logger.info(f"   Highly flexible residues (>{threshold:.2f} Å):")
        for res, flex in zip(high_flex_res, high_flex_rmsf):
            logger.info(f"     - Residue {res}: {flex:.3f} Å")



@time_tracker
def save_molecule_viewer_html(config):
    warnings.filterwarnings(
        "ignore", message="Found no information for attr: 'formalcharges'"
    )
    warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
    warnings.filterwarnings(
        "ignore", category=DeprecationWarning, module="MDAnalysis.coordinates.DCD"
    )

    emin_structure = config.get("path_emin_structure")
    md_trajectory = get_last_final_dcd(config)

    if not os.path.exists(emin_structure) or not os.path.exists(md_trajectory):
        logger.error(
            "❌ MD simulation files not found. Please run the quickrun() command first."
        )
    else:
        logger.info("Loading trajectory files...")

        u = mda.Universe(emin_structure, md_trajectory)
        system = u.select_atoms("not (resname HOH CL NA SOL WAT TIP3P)")

        logger.info(f"Creating animation with {len(u.trajectory)} frames...")

        temp_traj = "/tmp/trajectory_for_viz.pdb"
        system.write(temp_traj, frames="all")

        view = py3Dmol.view(width=800, height=600)
        with open(temp_traj) as f:
            view.addModelsAsFrames(f.read())

        # Set the viewer background to black
        view.setBackgroundColor("black")

        view.setStyle({"model": -1}, {"cartoon": {"color": "spectrum"}})
        view.setStyle(
            {"model": -1, "resn": ["UNK", "EPE"]},
            {"stick": {"radius": 0.3, "color": "red"}},
        )
        view.animate({"loop": "forward", "interval": 75})
        view.zoomTo()
        logger.info("\n=== 3D Trajectory Animation ===")
        view.show()

        html = view._make_html()

        # Create a more robust HTML modification with black background, centering, and title
        custom_html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>MD Trajectory Animation</title>
    <style>
        html, body {{
            background-color: black !important;
            margin: 0;
            padding: 0;
            font-family: Arial, sans-serif;
            color: white;
        }}
        .container {{
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            min-height: 100vh;
            padding: 20px;
        }}
        .title {{
            font-size: 24px;
            font-weight: bold;
            margin-bottom: 20px;
            text-align: center;
            color: white;
        }}
        .viewer-container {{
            border: 2px solid #333;
            border-radius: 10px;
            padding: 10px;
            background-color: black;
        }}
        .info {{
            margin-top: 15px;
            font-size: 14px;
            color: #ccc;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="title">MD Trajectory Animation</div>
        <div class="viewer-container">
            {html}
        </div>
        <div class="info">
            Total frames: {len(u.trajectory)} | Animation shows molecular dynamics simulation
        </div>
    </div>
</body>
</html>
"""

        with open(config.get("path_molviewer_html"), "w") as f:
            f.write(custom_html)

        logger.info(f"Location: {config.get('path_molviewer_html')}")
        logger.info(f"The HTML file contains the full animated trajectory")
        logger.info(f"Total frames in animation: {len(u.trajectory)}")


# --------------------------------------------------------------------------
# Main Function
# --------------------------------------------------------------------------


@time_tracker
def save_rmsd_rmsf_html(config):
    """Run RMSD and RMSF analysis on MD trajectory."""
    logger.info("========================================================")
    logger.info("📊 MD Trajectory Analysis")
    logger.info("========================================================")
    
    warnings.filterwarnings('ignore')
    
    u = None  # Initialize u to None for proper cleanup
    try:  
        # Load trajectory
        logger.info("Loading trajectory files...")
        last_dcd = get_last_final_dcd(config)
        u = mda.Universe(config.get("path_emin_structure"), last_dcd)
        
        # Calculate RMSD
        logger.info("Calculating RMSD...")
        rmsd_analysis, rmsd_ligand = _calculate_rmsd(u)
        
        # Calculate RMSF
        logger.info("Calculating RMSF...")
        rmsf_values, residue_ids = _calculate_rmsf(u)
        
        # Create interactive plot
        logger.info("Generating enhanced interactive plots...")
        fig = _create_interactive_plot(rmsd_analysis, rmsd_ligand, rmsf_values, residue_ids)
        
        # Save plot
        html_plot_path = config.get("path_rmsd_rmsf_html")
        pio.write_html(fig, file=html_plot_path, auto_open=False)
        logger.info(f"Interactive plot saved to: {html_plot_path}")
        
        # Create and save analysis summary table
        logger.info("Creating analysis summary table...")
        summary_html = _create_analysis_summary_table(u, rmsd_analysis, rmsd_ligand, rmsf_values, residue_ids)
        
        # Save summary table
        summary_path = html_plot_path.replace('.html', '_summary.html')
        with open(summary_path, 'w') as f:
            f.write(summary_html)
        logger.info(f"Analysis summary table saved to: {summary_path}")
        
        # Save analysis data
        _save_analysis_data(config, rmsd_analysis, rmsf_values, residue_ids, rmsd_ligand)
        
        # Log summary
        _log_analysis_summary(u, rmsd_analysis, rmsd_ligand, rmsf_values, residue_ids)
        
        logger.info("✅ Analysis Complete!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise
    finally:
        # Properly close the universe and trajectory to prevent DCDReader cleanup issues
        if u is not None:
            try:
                if hasattr(u, 'trajectory') and hasattr(u.trajectory, 'close'):
                    u.trajectory.close()
                # Clean up any remaining file handles
                if hasattr(u, '_trajectory'):
                    u._trajectory = None
            except (AttributeError, Exception) as cleanup_error:
                logger.warning(f"Warning during cleanup: {cleanup_error}")
                # Suppress cleanup errors to prevent masking the main error
                pass


# --------------------------------------------------------------------------
# Main Execution
# --------------------------------------------------------------------------


if __name__ == "__main__":
    config_filepath = "/Users/ingrid/Projects/EasyMD/easy-md/example/config/simulation_config.yaml"
    if config_filepath is not None:
        with open(config_filepath, "r") as f:
            config = yaml.safe_load(f)
    else:
        logger.error("No config file provided. Please provide a config file.")
        raise ValueError("No config file provided. Please provide a config file.")
    save_rmsd_rmsf_html(config)