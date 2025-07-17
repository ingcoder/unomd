import pytest
import os
import yaml
from pathlib import Path

# Test if all required packages are installed
def test_imports():
    """Test that all required packages can be imported."""
    try:
        import openmm
        import openff
        import numpy
        import mdtraj
    except ImportError as e:
        pytest.fail(f"Failed to import required package: {str(e)}")

def test_config_creation(tmp_path):
    """Test that configuration file can be created with basic settings."""
    from unomd.utils.config import create_config
    
    # Create a mock protein file
    protein_file = tmp_path / "test.pdb"
    protein_file.write_text("MOCK PDB CONTENT")
    
    # Create configuration
    config = create_config(
        protein_file=str(protein_file),
        project_dir=str(tmp_path),
        output_dir=str(tmp_path / "output")
    )
    
    # Check if config contains essential keys
    essential_keys = [
        'path_protein',
        'path_openmm_system',
        'path_openmm_topology',
        'platform_name',
        'platform_precision'
    ]
    
    for key in essential_keys:
        assert key in config, f"Missing essential key in config: {key}"
    
    # Check if output directory was created
    assert os.path.exists(tmp_path / "output"), "Output directory was not created"
    
    # Check if config file was saved
    config_file = tmp_path / "config" / "simulation_config.yaml"
    assert config_file.exists(), "Config file was not saved"
    
    # Verify config file content
    with open(config_file) as f:
        saved_config = yaml.safe_load(f)
    assert isinstance(saved_config, dict), "Saved config is not a valid YAML dictionary"

def test_directory_structure(tmp_path):
    """Test that the expected directory structure is created."""
    from unomd.utils.config import create_config
    
    # Create basic configuration
    protein_file = tmp_path / "test.pdb"
    protein_file.write_text("MOCK PDB CONTENT")
    
    config = create_config(
        protein_file=str(protein_file),
        project_dir=str(tmp_path)
    )
    
    # Check essential directories
    essential_dirs = [
        tmp_path / "output",
        tmp_path / "config"
    ]
    
    for directory in essential_dirs:
        assert directory.exists(), f"Directory not created: {directory}"
        assert directory.is_dir(), f"Path exists but is not a directory: {directory}"

def test_file_paths(tmp_path):
    """Test that file paths are correctly configured."""
    from unomd.utils.config import create_config
    
    # Create mock files
    protein_file = tmp_path / "test.pdb"
    protein_file.write_text("MOCK PDB CONTENT")
    
    config = create_config(
        protein_file=str(protein_file),
        project_dir=str(tmp_path)
    )
    
    # Check that paths are properly set
    assert config['path_protein'] == str(protein_file)
    assert config['path_base'] == str(tmp_path)
    assert config['path_openmm_system'].startswith(str(tmp_path))
    assert config['path_openmm_topology'].startswith(str(tmp_path))
    
    # Check that paths use correct separators for the OS
    for key, path in config.items():
        if key.startswith('path_') and isinstance(path, str) and path:  # Only check non-empty paths
            assert '\\\\' not in path, f"Invalid path separator in {key}: {path}"
            assert os.path.sep in path or path == str(protein_file), f"Missing path separator in {key}: {path}" 