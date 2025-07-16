# EasyMD Docker Setup

This guide explains how to run EasyMD using Docker containers for easy deployment and reproducible molecular dynamics simulations.

## Prerequisites

- Docker installed on your system
- Docker Compose (usually comes with Docker Desktop)
- For GPU support: NVIDIA Docker runtime

## Important Note About Dependencies

Due to conda channel availability issues, the Docker image is built with core dependencies only. **You'll need to manually install `ambertools` and `openff-toolkit` after building the container**. See the [Manual Dependency Installation](#manual-dependency-installation) section below.

## Quick Start

### 1. Build the Docker Image

```bash
# Build the image
docker-compose build easy-md

# Or build manually
docker build -t easy-md .
```

### 2. Install Missing Dependencies

The base image doesn't include `ambertools` and `openff-toolkit`. Install them manually:

```bash
# Start container
docker-compose up -d easy-md

# Install missing dependencies
docker-compose exec easy-md bash -c "
  conda activate easymd &&
  conda install -c conda-forge -c omnia ambertools &&
  conda install -c conda-forge openff-toolkit
"
```

### 3. Run a Simple Simulation

```bash
# Run a simulation with the example files
docker-compose exec easy-md easy-md /app/inputs/4W52.pdb --ligand_file /app/inputs/4w52_C_EPE.sdf --nsteps 1000
```

### 4. Check Results

```bash
# Check output files
docker-compose exec easy-md ls -la /app/outputs

# View logs
docker-compose exec easy-md ls -la /app/logs
```

## Manual Dependency Installation

### Option 1: Using conda-forge (Recommended)

```bash
# Start interactive session
docker-compose exec easy-md bash

# Activate environment
conda activate easymd

# Try installing from conda-forge
conda install -c conda-forge ambertools openff-toolkit
```

### Option 2: Using official channels

```bash
# For AmberTools (official channel)
conda install -c dacase ambertools-dac=25

# For OpenFF Toolkit
conda install -c conda-forge openff-toolkit
```

### Option 3: Manual installation from source

If conda installation fails, you can install from source:

```bash
# Inside the container
conda activate easymd

# Install AmberTools from source (requires more setup)
# See: https://ambermd.org/GetAmber.php

# Install OpenFF Toolkit via pip
pip install openff-toolkit
```

## Usage Options

### Option 1: Using Docker Compose (Recommended)

#### CPU-only simulation:

```bash
# Start container
docker-compose up -d easy-md

# Install dependencies (first time only)
docker-compose exec easy-md bash -c "
  conda activate easymd &&
  conda install -c conda-forge -c omnia ambertools openff-toolkit
"

# Run simulation
docker-compose exec easy-md easy-md /app/inputs/4W52.pdb --ligand_file /app/inputs/4w52_C_EPE.sdf --nsteps 5000 --platform CPU
```

#### GPU simulation (requires NVIDIA Docker):

```bash
# Start GPU container
docker-compose up -d easy-md-gpu

# Install dependencies (first time only)
docker-compose exec easy-md-gpu bash -c "
  conda activate easymd &&
  conda install -c conda-forge -c omnia ambertools openff-toolkit
"

# Run simulation with CUDA
docker-compose exec easy-md-gpu easy-md /app/inputs/4W52.pdb --ligand_file /app/inputs/4w52_C_EPE.sdf --nsteps 5000 --platform CUDA
```

#### Jupyter Notebook interface:

```bash
# Start Jupyter service
docker-compose up -d easy-md-jupyter

# Install dependencies (first time only)
docker-compose exec easy-md-jupyter bash -c "
  conda activate easymd &&
  conda install -c conda-forge -c omnia ambertools openff-toolkit
"

# Access at http://localhost:8888
```

### Option 2: Direct Docker Commands

```bash
# Build image
docker build -t easy-md .

# Run with volume mounts
docker run -v $(pwd)/example:/app/inputs:ro -v $(pwd)/outputs:/app/outputs -v $(pwd)/logs:/app/logs easy-md bash -c "
  conda activate easymd &&
  conda install -c conda-forge -c omnia ambertools openff-toolkit -y &&
  easy-md /app/inputs/4W52.pdb --ligand_file /app/inputs/4w52_C_EPE.sdf --nsteps 1000
"
```

## Creating a Custom Image with All Dependencies

If you want to create a custom image with all dependencies pre-installed:

```bash
# Start a container
docker run -it --name easy-md-setup easy-md bash

# Inside the container, install dependencies
conda activate easymd
conda install -c conda-forge -c omnia ambertools openff-toolkit -y

# Exit container
exit

# Commit the changes to a new image
docker commit easy-md-setup easy-md-complete

# Clean up
docker rm easy-md-setup

# Now you can use the complete image
docker run --rm easy-md-complete easy-md --help
```

## Command Line Options

The containerized `easy-md` command supports these options:

```bash
easy-md <protein_file> [options]

Required:
  protein_file          Path to protein PDB file

Optional:
  --ligand_file         Path to ligand file (SDF, MOL2, etc.)
  --nsteps              Number of MD steps (default: 1000)
  --platform            Platform: CPU or CUDA (default: CPU)
  --temperature         Temperature in Kelvin (default: 300.0)
  --output_dir          Output directory (default: /app/outputs)
```

## Directory Structure

The Docker container is organized as follows:

```
/app/
├── inputs/          # Input files (mounted from ./example)
├── outputs/         # Output files (mounted from ./outputs)
├── logs/           # Log files (mounted from ./logs)
├── notebooks/      # Jupyter notebooks (mounted from ./notebooks)
└── src/            # EasyMD source code
```

## Example Workflows

### Basic Protein-Ligand Simulation

```bash
# 1. Start container
docker-compose up -d easy-md

# 2. Install dependencies (first time only)
docker-compose exec easy-md bash -c "
  conda activate easymd &&
  conda install -c conda-forge -c omnia ambertools openff-toolkit -y
"

# 3. Run simulation
docker-compose exec easy-md easy-md \
  /app/inputs/4W52.pdb \
  --ligand_file /app/inputs/4w52_C_EPE.sdf \
  --nsteps 10000 \
  --temperature 310.0 \
  --platform CPU

# 4. Check results
docker-compose exec easy-md ls -la /app/outputs
```

### Using Custom Input Files

```bash
# 1. Place your files in the example directory
cp your_protein.pdb ./example/
cp your_ligand.sdf ./example/

# 2. Run simulation
docker-compose exec easy-md easy-md \
  /app/inputs/your_protein.pdb \
  --ligand_file /app/inputs/your_ligand.sdf \
  --nsteps 50000
```

## Troubleshooting

### Common Issues

1. **Missing dependencies error**: Make sure to install `ambertools` and `openff-toolkit` manually

   ```bash
   docker-compose exec easy-md bash -c "
     conda activate easymd &&
     conda install -c conda-forge -c omnia ambertools openff-toolkit -y
   "
   ```

2. **Permission errors**: Make sure the output directories are writable

   ```bash
   mkdir -p outputs logs
   chmod 777 outputs logs
   ```

3. **AmberTools installation fails**: Try different channels

   ```bash
   # Try official channel
   conda install -c dacase ambertools-dac=25

   # Or try omnia channel
   conda install -c omnia ambertools
   ```

4. **GPU not detected**: Ensure NVIDIA Docker runtime is installed
   ```bash
   # Check if nvidia-docker is available
   docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
   ```

### Debugging

```bash
# Check container logs
docker-compose logs easy-md

# Interactive debugging session
docker-compose exec easy-md bash

# Check conda environment
docker-compose exec easy-md conda info --envs
docker-compose exec easy-md conda list -n easymd

# Test import
docker-compose exec easy-md python -c "
  from easy_md.main.quickrun import quickrun
  print('✅ EasyMD imports successfully')
"
```

## Why Manual Installation is Needed

The `ambertools` package has complex dependencies and is not consistently available across conda channels in containerized environments. The channels that provide `ambertools` (`omnia`, `dacase`) often have:

1. **Inconsistent availability** in Docker builds
2. **Platform-specific issues** (especially for ARM64/Apple Silicon)
3. **Version conflicts** with other packages
4. **License restrictions** that make automated installation difficult

By providing a base image with core dependencies and manual installation instructions, you get:

- **Faster builds** for the base image
- **Flexibility** to choose the ambertools version that works for your system
- **Reliability** in different environments

## Performance Tips

1. **Use GPU** for large simulations when available
2. **Adjust memory** limits in Docker Desktop for better performance
3. **Use SSD storage** for better I/O performance
4. **Monitor resources** with `docker stats`
5. **Create custom images** with pre-installed dependencies for repeated use

## Cleanup

```bash
# Stop all services
docker-compose down

# Remove all containers and volumes
docker-compose down -v

# Remove images
docker rmi easy-md
```

This Docker setup provides a solid foundation for running EasyMD molecular dynamics simulations with manual dependency installation as needed.
