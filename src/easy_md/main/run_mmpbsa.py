"""Script for converting OpenMM/OpenFF molecular systems to AMBER format for MM/PBSA calculations.

This script provides functionality to:
1. Convert OpenMM/OpenFF molecular systems to AMBER format (prmtop/inpcrd)
2. Validate the resulting AMBER topology files
3. Handle force field parameters and bond types
4. Support both protein-ligand complexes and standalone molecules

The script includes robust error checking for:
- Bonds with missing type information
- Ligand parameter issues
- Counts number of molecules in the system

The conversion process uses OpenFF's Interchange to handle force field parameters
and ParmEd for the final AMBER format conversion. It supports the ff14SB force field
with OpenFF impropers for proteins and SMIRNOFF for ligands.

Note: The script is specifically designed to handle systems that will be used
for MM/PBSA calculations, ensuring proper handling of constrained bonds and
force field parameters.
"""

# Standard library imports
import logging
import os
from pathlib import Path
from typing import List, Set, Union
import subprocess
import yaml
import glob

# Custom imports
from easy_md.utils.fileparser import time_tracker
from easy_md.utils import amber_files
from easy_md.utils import info_logger
import logging

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------


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
    
def clean_mmpbsa_files():
    """Cleaning the mmpbsa files."""
    script_dir = Path(__file__).resolve().parent
    for file in script_dir.glob("_MMPBSA*"):
        try:
            file.unlink()
            logger.info(f"Deleted: {file}")
        except Exception as e:
            logger.error(f"Could not delete {file}: {e}")

def create_mmpbsa_input_file(config):
    """Creating the mmpbsa input file."""
    try:
        with open(config["path_mmpbsa_in"], "w") as f:
            f.write("&general\n")
            f.write("interval   = 10,\n")
            f.write("verbose    = 2,\n")
            f.write("endframe=100,\n")
            f.write("keep_files=2,\n")
            f.write("strip_mask= :HOH:CL:NA,\n")
            f.write("\n")
            f.write("/\n")
            f.write("&gb\n")
            f.write("igb=2, saltcon=0.150,\n")
            f.write("/\n")
            f.write("&pb\n")
            f.write("istrng=0.150,inp=2, radiopt=0, prbrad=1.4,\n")
            f.write("/\n")
        logger.info(f"✅ Successfully created mmpbsa input file: {config['path_mmpbsa_in']}")
    except Exception as e:
        logger.error(f"Error creating mmpbsa input file: {e}")
        raise 


# ------------------------------------------------------------------------------
# Amber Tools
# ------------------------------------------------------------------------------


@time_tracker
def run_ante_mmpbsa(config) -> None:
    """Preparing the amber input topology files for the MMPBSA calculation."""

    command = [
        "ante-MMPBSA.py",
        "-p",
        config["path_amber_solvated"],
        "-c",
        config["path_amber_complex"],
        "-r",
        config["path_amber_receptor"],
        "-l",
        config["path_amber_ligand"],
        "-s",
        ":HOH:CL:NA",
        "-n",
        ":UNK",
        "--radii",
        "mbondi2",
    ]

    logger.info("Running ante-MMPBSA.py with the following command:")
    logger.info(" ".join(command))

    try:
        subprocess.run(command, check=True)
        logger.info("✅ ante-MMPBSA.py ran successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ ante-MMPBSA.py failed with return code {e.returncode}")
    except FileNotFoundError:
        logger.error("❌ ante-MMPBSA.py not found. Is it installed and in your PATH?")

@time_tracker
def run_mmpbsa(config) -> None:
    """Running the MMPBSA calculation."""

    last_final_dcd = get_last_final_dcd(config)

    command = [
        "MMPBSA.py",
        "-O",
        "-i",
        config["path_mmpbsa_in"],
        "-o",
        config["path_mmpbsa_results"],
        "-sp",
        config["path_amber_solvated"],
        "-cp",
        config["path_amber_complex"],
        "-rp",
        config["path_amber_receptor"],
        "-lp",
        config["path_amber_ligand"],
        "-y",
        last_final_dcd,
    ]

    logger.info("Running mmpbsa.py with the following command:")
    logger.info(" ".join(command))

    try:
        subprocess.run(command, check=True)
        logger.info("✅ mmpbsa.py ran successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ mmpbsa.py failed with return code {e.returncode}")
    except FileNotFoundError:
        logger.error("❌ mmpbsa.py not found. Is it installed and in your PATH?")

# ------------------------------------------------------------------------------
# Main Function
# ------------------------------------------------------------------------------


def main(config_filepath) -> None:
    logger.info("========================================================")
    logger.info("🛠️ MM/PBSA Calculation")
    logger.info("========================================================")

    if config_filepath is not None:
        with open(config_filepath, "r") as f:
            config = yaml.safe_load(f)
    else:
        logger.error("No config file provided. Please provide a config file.")
        raise ValueError("No config file provided. Please provide a config file.")

    create_mmpbsa_input_file(config)

    try:
        amber_files.export_prmtop(config)
        run_ante_mmpbsa(config)
        run_mmpbsa(config)
        clean_mmpbsa_files(config)
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        raise
    except PermissionError as e:
        logger.error(f"Permission error: {e}")
        raise
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


if __name__ == "__main__":
    config_filepath = "/Users/ingrid/Projects/EasyMD/easy-md/example/config/simulation_config.yaml"
    main(config_filepath)
