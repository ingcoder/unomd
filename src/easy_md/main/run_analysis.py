"""
RMSD and RMSF analysis of MD trajectories.

This module provides functionality to analyze molecular dynamics trajectories
by calculating Root Mean Square Deviation (RMSD) and Root Mean Square Fluctuation (RMSF)
with interactive visualization using Plotly.
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
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# Custom imports
from easy_md.utils.fileparser import time_tracker
from easy_md.utils import info_logger
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
    """Create interactive Plotly visualization."""
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Root Mean Square Deviation (RMSD)", "Root Mean Square Fluctuation (RMSF)")
    )
    
    # RMSD plot
    time = rmsd_analysis.rmsd[:, 1]
    rmsd_backbone = rmsd_analysis.rmsd[:, 2]
    
    fig.add_trace(
        go.Scatter(
            x=time, y=rmsd_backbone, mode='lines', name='Protein Backbone',
            line=dict(color='blue', width=2),
            hovertemplate="Time: %{x:.2f} ps<br>RMSD: %{y:.2f} Å<extra></extra>"
        ),
        row=1, col=1
    )
    
    # Add ligand RMSD if available
    if rmsd_ligand is not None:
        time_ligand = rmsd_ligand.rmsd[:, 1]
        rmsd_ligand_values = rmsd_ligand.rmsd[:, 2]
        fig.add_trace(
            go.Scatter(
                x=time_ligand, y=rmsd_ligand_values, mode='lines', name='Ligand',
                line=dict(color='red', width=2),
                hovertemplate="Time: %{x:.2f} ps<br>RMSD: %{y:.2f} Å<extra></extra>"
            ),
            row=1, col=1
        )
    
    # RMSF plot
    fig.add_trace(
        go.Scatter(
            x=residue_ids, y=rmsf_values, mode='lines', name='RMSF',
            line=dict(color='green', width=2), fill='tozeroy', fillcolor='rgba(0,255,0,0.2)',
            hovertemplate="Residue: %{x}<br>RMSF: %{y:.2f} Å<extra></extra>"
        ),
        row=2, col=1
    )
    
    # Add flexibility threshold
    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    threshold = mean_rmsf + 1.5 * std_rmsf
    
    fig.add_hline(
        y=threshold, line_dash="dash", line_color="red",
        annotation_text=f'Flexibility threshold ({threshold:.2f} Å)',
        annotation_position="top right", row=2, col=1
    )
    
    # Highlight highly flexible residues
    high_flex_mask = rmsf_values > threshold
    if np.any(high_flex_mask):
        high_flex_res = residue_ids[high_flex_mask]
        high_flex_rmsf = rmsf_values[high_flex_mask]
        fig.add_trace(
            go.Scatter(
                x=high_flex_res, y=high_flex_rmsf, mode='markers',
                name='Highly Flexible Residues',
                marker=dict(color='red', size=8, symbol='circle'),
                hovertemplate="Residue: %{x}<br>RMSF: %{y:.2f} Å<extra></extra>"
            ),
            row=2, col=1
        )
    
    # Update layout
    fig.update_xaxes(title_text='Time (ps)', row=1, col=1)
    fig.update_yaxes(title_text='RMSD (Å)', row=1, col=1)
    fig.update_xaxes(title_text='Residue Number', row=2, col=1)
    fig.update_yaxes(title_text='RMSF (Å)', row=2, col=1)
    
    fig.update_layout(
        title_text='MD Trajectory Analysis: RMSD and RMSF',
        title_x=0.5, height=800, showlegend=True, hovermode="x unified"
    )
    
    return fig


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


# --------------------------------------------------------------------------
# Main Function
# --------------------------------------------------------------------------


@time_tracker
def main(config):
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
        logger.info("Generating interactive plots...")
        fig = _create_interactive_plot(rmsd_analysis, rmsd_ligand, rmsf_values, residue_ids)
        
        # Save plot
        html_plot_path = config.get("path_rmsd_rmsf_html")
        pio.write_html(fig, file=html_plot_path, auto_open=False)
        logger.info(f"Interactive plot saved to: {html_plot_path}")
        
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
    main(config)
