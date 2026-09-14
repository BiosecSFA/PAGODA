#!/bin/bash

# Setup script for Overlapping Gene Design Pipeline - HPC/Linux with CUDA
# This script creates a conda environment and installs all dependencies

set -e  # Exit on error

echo "================================================================================"
echo "Overlapping Gene Design Pipeline - HPC/Linux Setup"
echo "================================================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Step 1: Check if conda is installed
echo "Step 1: Checking for conda/mamba..."
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "✓ Found mamba: $(mamba --version | head -1)"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo "✓ Found conda: $(conda --version)"
else
    echo "✗ Neither conda nor mamba found!"
    echo ""
    echo "Please install Miniconda or Mambaforge:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    echo "  https://github.com/conda-forge/miniforge"
    echo ""
    echo "For Linux:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh"
    echo ""
    exit 1
fi
echo ""

# Step 2: Check for CUDA module
echo "Step 2: Checking for CUDA..."
if command -v module &> /dev/null; then
    echo "  Module system detected"
    if module avail cuda 2>&1 | grep -q cuda; then
        echo "  ✓ CUDA modules available"
        module load cuda 
    else
        echo "  ⚠ CUDA modules not found"
        echo "  CCMpred compilation may fail without CUDA"
    fi
else
    echo "  ⚠ Module system not detected"
    echo "  Ensure CUDA toolkit is installed for CCMpred"
fi
echo ""

# Step 3: Create conda environment
echo "Step 3: Creating conda environment in conda/ directory..."
CONDA_ENV_DIR="$SCRIPT_DIR/conda"
if [ -d "$CONDA_ENV_DIR/bin" ]; then
    echo "  Environment already exists at conda/"
    read -p "  Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$CONDA_ENV_DIR"
        ${CONDA_CMD} env create -f environment_hpc.yml --prefix "$CONDA_ENV_DIR"
    fi
else
    ${CONDA_CMD} env create -f environment_hpc.yml --prefix "$CONDA_ENV_DIR"
fi
echo "✓ Conda environment created at conda/"
echo ""

# Step 4: Activate environment
echo "Step 4: Activating conda environment..."
# Initialize conda for this shell session
eval "$(conda shell.bash hook 2>/dev/null)" || eval "$(conda shell.zsh hook 2>/dev/null)" || true
source activate "$CONDA_ENV_DIR"
echo "✓ Environment activated"
echo ""

# Step 5: Install HMMER 3.4 from source
echo "Step 5: Installing HMMER 3.4..."
HMMER_DIR="$CONDA_ENV_DIR/external/hmmer"
if [ -f "$HMMER_DIR/bin/phmmer" ]; then
    echo "  ✓ HMMER already installed at conda/external/hmmer"
else
    echo "  Building HMMER from source..."
    mkdir -p "$CONDA_ENV_DIR/external"
    cd "$CONDA_ENV_DIR/external"

    # Download HMMER 3.4
    if [ ! -f "hmmer-3.4.tar.gz" ]; then
        echo "    Downloading HMMER 3.4..."
        wget -q http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz || \
            curl -L -o hmmer-3.4.tar.gz http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
    fi

    # Extract
    if [ ! -d "hmmer-3.4" ]; then
        echo "    Extracting..."
        tar -xzf hmmer-3.4.tar.gz
    fi

    # Build and install
    cd hmmer-3.4
    echo "    Configuring..."
    ./configure --prefix="$HMMER_DIR" > /dev/null 2>&1
    echo "    Building (this may take a few minutes)..."
    make -j4 > /dev/null 2>&1
    echo "    Installing..."
    make install > /dev/null 2>&1

    # Build Easel tools
    cd easel
    echo "    Building Easel tools..."
    make install > /dev/null 2>&1

    cd "$SCRIPT_DIR"
    echo "  ✓ HMMER 3.4 installed to conda/external/hmmer"
fi
echo ""

# Step 6: Install CCMpred with CUDA
echo "Step 6: Installing CCMpred (requires CUDA)..."
CCMPRED_DIR="$CONDA_ENV_DIR/external/ccmpred"
if [ -f "$CCMPRED_DIR/bin/ccmpred" ]; then
    echo "  ✓ CCMpred already installed at conda/external/ccmpred"
else
    echo "  Building CCMpred from source..."
    mkdir -p "$CONDA_ENV_DIR/external"
    cd "$CONDA_ENV_DIR/external"

    # Clone CCMpred
    if [ ! -d "CCMpred" ]; then
        echo "    Cloning CCMpred..."
        git clone --recursive https://github.com/soedinglab/CCMpred.git > /dev/null 2>&1
    fi

    cd CCMpred

    # Patch CMakeLists.txt for newer CMake versions
    echo "    Patching for CMake compatibility..."
    if [ -f "lib/libconjugrad/CMakeLists.txt" ]; then
        # Remove the problematic cmake_policy line
        sed -i.bak '/cmake_policy(SET CMP0022 OLD)/d' lib/libconjugrad/CMakeLists.txt
        # Update cmake_minimum_required to 3.10
        sed -i.bak 's/cmake_minimum_required(VERSION [0-9.]*)/cmake_minimum_required(VERSION 3.10)/' lib/libconjugrad/CMakeLists.txt
    fi

    # Build with cmake
    echo "    Running cmake..."
    if cmake .; then
        echo "    Building (this may take several minutes)..."
        if make -j4; then
            # Install
            mkdir -p "$CCMPRED_DIR/bin"

            if [ -f "bin/ccmpred" ]; then
                cp bin/ccmpred "$CCMPRED_DIR/bin/"
                echo "  ✓ CCMpred installed to conda/external/ccmpred/bin"
            elif [ -f "ccmpred" ]; then
                cp ccmpred "$CCMPRED_DIR/bin/"
                echo "  ✓ CCMpred installed to conda/external/ccmpred/bin"
            else
                echo "  ⚠ CCMpred built but executable not found"
                echo "     Check build directory for errors"
            fi
        else
            echo "  ✗ CCMpred compilation failed"
            echo "     Check error messages above"
            echo "     Common issues:"
            echo "       - CUDA not loaded: module load cuda"
            echo "       - Wrong CUDA version"
            echo "       - Missing cmake or build tools"
        fi
    else
        echo "  ✗ CMake configuration failed"
        echo "     Check error messages above"
        echo "     Common issues:"
        echo "       - CUDA not loaded: module load cuda"
        echo "       - Missing cmake: conda install cmake"
    fi

    cd "$SCRIPT_DIR"
fi
echo ""

# Step 7: Clone pycameox
echo "Step 7: Installing pycameox..."
PYCAMEOX_DIR="$CONDA_ENV_DIR/external/pycameox"
if [ -d "$PYCAMEOX_DIR" ]; then
    echo "  ✓ pycameox already exists at conda/external/pycameox"
else
    mkdir -p "$CONDA_ENV_DIR/external"
    cd "$CONDA_ENV_DIR/external"
    echo "  Cloning pycameox from GitHub..."
    git clone https://github.com/BiosecSFA/pycameox.git > /dev/null 2>&1
    cd pycameox
    echo "  Checking out commit aa0c040..."
    git checkout aa0c040 > /dev/null 2>&1
    cd "$SCRIPT_DIR"
    echo "  ✓ pycameox installed at conda/external/pycameox"
fi
echo ""

# Step 8: Clone potts
echo "Step 8: Installing potts..."
POTTS_DIR="$CONDA_ENV_DIR/external/potts"
if [ -d "$POTTS_DIR" ]; then
    echo "  ✓ potts already exists at conda/external/potts"
else
    mkdir -p "$CONDA_ENV_DIR/external"
    cd "$CONDA_ENV_DIR/external"
    echo "  Cloning potts from GitHub..."
    git clone https://github.com/hnisonoff/potts.git > /dev/null 2>&1
    cd potts
    echo "  Checking out commit 1b32512..."
    git checkout 1b32512 > /dev/null 2>&1
    cd "$SCRIPT_DIR"
    echo "  ✓ potts installed at conda/external/potts"
fi
echo ""

# Step 9: Create activation script
echo "Step 9: Creating activation script..."
cat > activate_release.sh << 'ACTIVATE_EOF'
#!/bin/bash
# Activation script for Overlapping Gene Design Pipeline Release

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Activate conda environment from local conda/ directory
conda activate "$SCRIPT_DIR/conda"

# CRITICAL: Disable user site-packages to prevent conflicts
# This prevents Python from loading packages from ~/.local/lib/python*/site-packages
export PYTHONNOUSERSITE=1

# Add HMMER and Easel to PATH (from conda/external/)
if [ -d "$SCRIPT_DIR/conda/external/hmmer/bin" ]; then
    export PATH="$SCRIPT_DIR/conda/external/hmmer/bin:$PATH"
fi

# Add ccmpred to PATH (from conda/external/)
if [ -f "$SCRIPT_DIR/conda/external/ccmpred/bin/ccmpred" ]; then
    export PATH="$SCRIPT_DIR/conda/external/ccmpred/bin:$PATH"
fi

# Add repository to PYTHONPATH
export PYTHONPATH="$SCRIPT_DIR:$PYTHONPATH"

# Add pycameox to PYTHONPATH (from conda/external/)
export PYTHONPATH="$SCRIPT_DIR/conda/external/pycameox/src:$PYTHONPATH"

# Add potts to PYTHONPATH (from conda/external/)
if [ -d "$SCRIPT_DIR/conda/external/potts" ]; then
    export PYTHONPATH="$SCRIPT_DIR/conda/external/potts:$PYTHONPATH"
fi

echo "✓ Overlap Design environment activated"
echo ""
echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""
if [ -f "$SCRIPT_DIR/conda/external/hmmer/bin/phmmer" ]; then
    echo "HMMER: $($SCRIPT_DIR/conda/external/hmmer/bin/phmmer -h 2>&1 | grep "HMMER" | head -n1)"
else
    echo "HMMER: Not found"
fi
if command -v ccmpred > /dev/null; then
    echo "CCMpred: $(which ccmpred)"
else
    echo "CCMpred: Not found (load CUDA module and rebuild if needed)"
fi
echo ""
echo "Example commands:"
echo ""
echo "  # Design overlaps with controlled overlap range (10-30 aa)"
echo "  python scripts/dpalign.py purE hisI 1 --min-overlap 10 --max-overlap 30"
echo ""
echo "  # Test new gene against conditionally essential genes"
echo "  python test_new_gene.py newgene.fasta \\"
echo "      --database /p/lustre1/xu26/datangle/dbs/uniref90.fasta.gz \\"
echo "      --cds-file data/split_cds/newgene.fasta \\"
echo "      --min-overlap 10 --max-overlap 30"
echo ""
echo "  # Visualize Potts model couplings"
echo "  python scripts/plot_ccmpred_circos.py data/raw/purE.raw purE_circos.png \\"
echo "      --fasta data/split_fasta/purE.fasta"
echo ""
echo "To run tests:"
echo "  bash test_pipeline.sh"
echo ""
ACTIVATE_EOF

chmod +x activate_release.sh
echo "  ✓ Created activate_release.sh"
echo ""

# Summary
echo "================================================================================"
echo "Setup Complete!"
echo "================================================================================"
echo ""
echo "Installed components:"
echo "  ✓ Conda environment at conda/"
echo "  ✓ Python $(python --version 2>&1 | awk '{print $2}')"
echo "  ✓ HMMER 3.4 at conda/external/hmmer/"
echo "  ✓ pycameox at conda/external/pycameox/"
echo "  ✓ potts at conda/external/potts/"
if [ -f "$CCMPRED_DIR/bin/ccmpred" ]; then
    echo "  ✓ CCMpred at conda/external/ccmpred/"
else
    echo "  ⚠ CCMpred: Build failed or CUDA not available"
fi
echo ""
echo "Example genes included:"
echo "  • purE (169 aa) - Small, experimental validation"
echo "  • hisI (202 aa) - Small-medium, experimental validation"
echo "  • ilvE (309 aa) - Medium, experimental validation"
echo "  • gltA (427 aa) - Medium-large, experimental validation"
echo ""
echo "Next steps:"
echo "  1. Load CUDA module (if not already loaded):"
echo "       module load cuda"
echo ""
echo "  2. Download UniRef90 database (required for new gene model generation):"
echo "       mkdir -p data/db"
echo "       cd data/db"
echo "       wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
echo "       cd ../.."
echo "       Note: This is a ~20GB download (50GB uncompressed)"
echo ""
echo "  3. Activate the environment:"
echo "       source activate_release.sh"
echo ""
echo "  4. Run tests:"
echo "       bash test_pipeline.sh"
echo ""
echo "  5. Try an example:"
echo "       python scripts/dpalign.py purE hisI 1 --min-overlap 10 --max-overlap 30"
echo ""
echo "================================================================================"
