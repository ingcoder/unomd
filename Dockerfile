# Use Ubuntu 22.04 as base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge for x86_64, which comes with mamba and conda-forge pre-configured
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

# Add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Copy environment and project files
COPY environment.yaml pyproject.toml ./
COPY src/ ./src/
COPY README.md LICENSE ./

# Create conda environment from environment.yaml using mamba
RUN mamba env create -f environment.yaml && mamba clean -afy

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "easymd", "/bin/bash", "-c"]

# Install ambertools and openff-toolkit
RUN mamba install -c conda-forge -c omnia ambertools openff-toolkit -y

# Install the package in development mode
RUN pip install -e .

# Set the default command to open an interactive bash shell
CMD ["/bin/bash"]