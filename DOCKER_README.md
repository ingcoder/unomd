# Docker Quick Start

## Prerequisites

- Docker installed on your system
- Docker Compose installed

## Build the Image

```bash
docker build --platform=linux/amd64 -t unomd .
```

## Run the Container

```bash
docker run -it unomd
```

## Basic Usage

Once inside the container, you can run UnoMD simulations:

```bash
# Example simulation
unomd protein.pdb --ligand_file ligand.sdf --nsteps 1000
```

## Output

Results will be saved to the `outputs/` directory, which is mounted to your host system.
