#!/bin/bash

# UnoMD Dependency Installation Script
# This script installs ambertools and openff-toolkit in the Docker container

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if running inside Docker container
if [ ! -f /.dockerenv ]; then
    print_error "This script should be run inside the Docker container"
    print_info "Run: docker-compose exec uno-md bash"
    print_info "Then: ./install_dependencies.sh"
    exit 1
fi

# Activate conda environment
print_info "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate unomd

# Function to try installing from different channels
install_ambertools() {
    print_info "Installing AmberTools..."
    
    # Try conda-forge first
    print_info "Trying conda-forge channel..."
    if conda install -c conda-forge ambertools -y; then
        print_info "✅ AmberTools installed successfully from conda-forge"
        return 0
    fi
    
    # Try omnia channel
    print_info "Trying omnia channel..."
    if conda install -c omnia ambertools -y; then
        print_info "✅ AmberTools installed successfully from omnia"
        return 0
    fi
    
    # Try official dacase channel
    print_info "Trying official dacase channel..."
    if conda install -c dacase ambertools-dac=25 -y; then
        print_info "✅ AmberTools installed successfully from dacase"
        return 0
    fi
    
    print_error "Failed to install AmberTools from all channels"
    return 1
}

install_openff_toolkit() {
    print_info "Installing OpenFF Toolkit..."
    
    # Try conda-forge
    if conda install -c conda-forge openff-toolkit -y; then
        print_info "✅ OpenFF Toolkit installed successfully"
        return 0
    fi
    
    # Try pip as fallback
    print_info "Trying pip installation..."
    if pip install openff-toolkit; then
        print_info "✅ OpenFF Toolkit installed successfully via pip"
        return 0
    fi
    
    print_error "Failed to install OpenFF Toolkit"
    return 1
}

# Test imports
test_imports() {
    print_info "Testing imports..."
    
    python -c "
import sys
try:
    from openff.toolkit import ForceField
    print('✅ OpenFF Toolkit import successful')
except ImportError as e:
    print('❌ OpenFF Toolkit import failed:', e)
    sys.exit(1)

try:
    from easy_md.main.quickrun import quickrun
    print('✅ EasyMD import successful')
except ImportError as e:
    print('❌ EasyMD import failed:', e)
    sys.exit(1)

print('🎉 All dependencies are working correctly!')
"
}

# Main installation process
main() {
    print_info "Starting EasyMD dependency installation..."
    
    # Install AmberTools
    if ! install_ambertools; then
        print_error "AmberTools installation failed"
        exit 1
    fi
    
    # Install OpenFF Toolkit
    if ! install_openff_toolkit; then
        print_error "OpenFF Toolkit installation failed"
        exit 1
    fi
    
    # Test imports
    test_imports
    
    print_info "🎉 All dependencies installed successfully!"
    print_info "You can now run EasyMD simulations."
}

# Run main function
main "$@" 