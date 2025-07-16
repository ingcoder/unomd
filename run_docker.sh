#!/bin/bash

# EasyMD Docker Helper Script
# This script helps you build and run EasyMD in Docker containers

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if Docker is installed
check_docker() {
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install Docker first."
        exit 1
    fi
    
    if ! command -v docker-compose &> /dev/null; then
        print_error "Docker Compose is not installed. Please install Docker Compose first."
        exit 1
    fi
}

# Create necessary directories
setup_directories() {
    print_info "Setting up directories..."
    mkdir -p outputs logs
    chmod 755 outputs logs
}

# Build Docker image
build() {
    print_info "Building EasyMD Docker image..."
    docker-compose build easy-md
    print_info "Build completed successfully!"
}

# Run a simple test
test_run() {
    print_info "Running a test simulation..."
    setup_directories
    
    # Start container
    docker-compose up -d easy-md
    
    # Run test simulation
    docker-compose exec easy-md easy-md /app/inputs/4W52.pdb --ligand_file /app/inputs/4w52_C_EPE.sdf --nsteps 100
    
    print_info "Test completed! Check the outputs directory for results."
}

# Start interactive session
interactive() {
    print_info "Starting interactive session..."
    setup_directories
    docker-compose up -d easy-md
    docker-compose exec easy-md bash
}

# Start Jupyter notebook
jupyter() {
    print_info "Starting Jupyter notebook..."
    setup_directories
    docker-compose up -d easy-md-jupyter
    print_info "Jupyter notebook is available at: http://localhost:8888"
}

# Start GPU container
gpu() {
    print_info "Starting GPU container..."
    setup_directories
    
    # Check if nvidia-docker is available
    if ! docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi &> /dev/null; then
        print_warning "NVIDIA Docker runtime not detected. GPU acceleration may not work."
    fi
    
    docker-compose up -d easy-md-gpu
    print_info "GPU container started. You can now run GPU-accelerated simulations."
}

# Clean up
cleanup() {
    print_info "Cleaning up Docker resources..."
    docker-compose down -v
    docker rmi easy-md 2>/dev/null || true
    print_info "Cleanup completed!"
}

# Show usage
usage() {
    cat << EOF
EasyMD Docker Helper Script

Usage: $0 [COMMAND]

Commands:
    build           Build the Docker image
    test            Run a test simulation
    interactive     Start interactive session
    jupyter         Start Jupyter notebook server
    gpu             Start GPU-enabled container
    cleanup         Clean up Docker resources
    help            Show this help message

Examples:
    $0 build                    # Build the image
    $0 test                     # Run a quick test
    $0 interactive              # Start interactive session
    $0 jupyter                  # Start Jupyter notebook
    $0 gpu                      # Start GPU container

For more detailed usage, see DOCKER_README.md
EOF
}

# Main script logic
main() {
    case "${1:-help}" in
        build)
            check_docker
            build
            ;;
        test)
            check_docker
            test_run
            ;;
        interactive)
            check_docker
            interactive
            ;;
        jupyter)
            check_docker
            jupyter
            ;;
        gpu)
            check_docker
            gpu
            ;;
        cleanup)
            check_docker
            cleanup
            ;;
        help|--help|-h)
            usage
            ;;
        *)
            print_error "Unknown command: $1"
            usage
            exit 1
            ;;
    esac
}

# Run main function
main "$@" 