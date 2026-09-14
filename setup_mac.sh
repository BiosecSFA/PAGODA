#!/bin/bash

# Setup script for Overlapping Gene Design Pipeline - macOS
# This script creates a conda environment and installs all dependencies

set -e  # Exit on error

echo "================================================================================"
echo "Overlapping Gene Design Pipeline - macOS Setup"
echo "================================================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Step 1: Check if conda is installed
echo "Step 1: Checking for conda..."
if ! command -v conda &> /dev/null; then
    echo "✗ Conda not found!"
    echo ""
    echo "Please install Miniconda or Anaconda:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    echo ""
    echo "For macOS:"
    echo "  curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    echo "  bash Miniconda3-latest-MacOSX-x86_64.sh"
    echo ""
    echo "For macOS arm64:"
    echo "  https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    echo "  bash Miniconda3-latest-MacOSX-arm64.sh"
    echo ""    
    exit 1
fi
echo "✓ Found conda: $(conda --version)"
echo ""

# Step 2: Create conda environment
echo "Step 2: Creating conda environment in conda/ directory..."
CONDA_ENV_DIR="$SCRIPT_DIR/conda"
if [ -d "$CONDA_ENV_DIR/bin" ]; then
    echo "  Environment already exists at conda/"
    read -p "  Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$CONDA_ENV_DIR"
        conda env create -f environment_mac.yml --prefix "$CONDA_ENV_DIR"
    fi
else
    conda env create -f environment_mac.yml --prefix "$CONDA_ENV_DIR"
fi
echo "✓ Conda environment created at conda/"
echo ""

# Step 3: Activate environment
echo "Step 3: Activating conda environment..."
# Initialize conda for this shell session
eval "$(conda shell.bash hook 2>/dev/null)" || eval "$(conda shell.zsh hook 2>/dev/null)" || true
source activate "$CONDA_ENV_DIR"
echo "✓ Environment activated"
echo ""

# Step 4: Install HMMER 3.4 from source
echo "Step 4: Installing HMMER 3.4..."
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

# Step 5: Install CCMpred (CPU-only mode for macOS)
echo "Step 5: Installing CCMpred (CPU-only, no CUDA)..."
CCMPRED_DIR="$CONDA_ENV_DIR/external/ccmpred"
if [ -f "$CCMPRED_DIR/bin/ccmpred" ]; then
    echo "  ✓ CCMpred already installed at conda/external/ccmpred"
else
    echo "  Building CCMpred from source (CPU-only mode)..."
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

    # Build with cmake (CUDA disabled for macOS)
    echo "    Running cmake with CUDA disabled..."
    if cmake -DWITH_CUDA=OFF .; then
        echo "    Building (this may take several minutes)..."
        if make -j4; then
            # Install
            mkdir -p "$CCMPRED_DIR/bin"

            if [ -f "bin/ccmpred" ]; then
                echo "  ✓ CCMpred (CPU-only) installed to conda/external/ccmpred/bin"
            elif [ -f "ccmpred" ]; then
                cp ccmpred "$CCMPRED_DIR/bin/"
                echo "  ✓ CCMpred (CPU-only) installed to conda/external/ccmpred/bin"
            else
                echo "  ⚠ CCMpred built but executable not found"
                echo "     Check build directory for details"
            fi
        else
            echo "  ✗ CCMpred compilation failed"
            echo "     Check error messages above"
            echo "     You can still use pre-computed Potts models"
        fi
    else
        echo "  ✗ CMake configuration failed"
        echo "     Check error messages above"
        echo "     You can still use pre-computed Potts models"
    fi

    cd "$SCRIPT_DIR"
fi
echo ""
echo "  Note: CCMpred on macOS runs in CPU-only mode (slower than GPU)"
echo "        For faster Potts model generation, use HPC with CUDA"
echo ""

# Step 6: Clone pycameox
echo "Step 6: Installing pycameox..."
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

# Step 7: Clone potts
echo "Step 7: Installing potts..."
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

# Step 8: Create activation script
echo "Step 8: Creating activation script..."
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

# Add CCMpred to PATH (from conda/external/)
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
if command -v ccmpred > /dev/null 2>&1; then
    echo "CCMpred: $(which ccmpred) (CPU-only mode)"
else
    echo "CCMpred: Not found (use pre-computed .raw files or rebuild)"
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
    echo "  ✓ CCMpred (CPU-only) at conda/external/ccmpred/"
else
    echo "  ⚠ CCMpred: Build failed (use pre-computed Potts models)"
fi
echo ""
echo "Example genes included:"
echo "  • purE (169 aa) - Small, experimental validation"
echo "  • hisI (202 aa) - Small-medium, experimental validation"
echo "  • ilvE (309 aa) - Medium, experimental validation"
echo "  • gltA (427 aa) - Medium-large, experimental validation"
echo ""
echo "Next steps:"
echo "  1. Download UniRef90 database (required for new gene model generation):"
echo "       mkdir -p data/db"
echo "       cd data/db"
echo "       curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
echo "       cd ../.."
echo "       Note: This is a ~20GB download (50GB uncompressed)"
echo "       Note: On macOS, CCMpred is not available, so Potts models cannot be generated"
echo ""
echo "  2. Activate the environment:"
echo "       source activate_release.sh"
echo ""
echo "  3. Run tests:"
echo "       bash test_pipeline.sh"
echo ""
echo "  4. Try an example:"
echo "       python scripts/dpalign.py purE hisI 1 --min-overlap 10 --max-overlap 30"
echo ""
echo "================================================================================"
