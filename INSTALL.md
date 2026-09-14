# Installation Guide

## Quick Installation

### macOS
```bash
bash setup_mac.sh
```

### HPC/Linux
```bash
bash setup_hpc.sh
```

## What Gets Installed

All components are installed within the `release/` directory:

```
release/
├── conda/                    # Self-contained conda environment (~2-3GB)
│   ├── bin/                 # Python 3.11, conda packages
│   ├── lib/
│   └── ...
├── external/
│   ├── hmmer/               # HMMER 3.4 (~50MB)
│   ├── ccmpred/             # CCMpred (HPC only, ~10MB)
│   ├── pycameox/            # pyCameoX package
│   └── potts/               # Potts package
└── data/                    # Example gene data (~2GB)
```

**Total Size**: ~5GB after installation

## Installation Steps Explained

### Step 1: Check for Conda

The setup script checks if conda or mamba is installed. If not, it will provide download instructions.

### Step 2: Create Conda Environment

Creates a local conda environment at `conda/` with:
- Python 3.11
- Scientific packages: biotite, biopython, numpy, pandas, pytorch, scipy
- Visualization: matplotlib, seaborn
- Build tools: cmake, make, git

**Time**: 5-10 minutes depending on internet speed

### Step 3: Build HMMER

Downloads and compiles HMMER 3.4 from source into `external/hmmer/`:
- phmmer: Search for homologs
- hmmbuild: Build profile HMMs
- esl-reformat: Format conversion

**Time**: 2-5 minutes

### Step 4: Build CCMpred (HPC only)

Compiles CCMpred with CUDA support for Potts model training.

**Time**: 3-5 minutes
**Requires**: CUDA toolkit (module load cuda)

### Step 5: Install pycameox

Clones pycameox from GitHub (commit aa0c040) into `external/pycameox/`.

### Step 6: Install potts

Clones potts package from GitHub (commit 1b32512) into `external/potts/`.

### Step 7: Create Activation Script

Generates `activate_release.sh` that:
- Activates the local conda environment
- Adds HMMER/CCMpred to PATH
- Sets PYTHONPATH for custom packages

## Verifying Installation

After setup completes:

```bash
# Activate environment
source activate_release.sh

# Check Python
python --version  # Should be 3.11.x

# Check tools
phmmer -h | head -1  # Should show HMMER 3.4
ccmpred -h           # (HPC only) Should show CCMpred help

# Run tests
bash test_pipeline.sh
```

## Troubleshooting

### Conda Not Found

**macOS:**
```bash
# Download Miniconda
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Or with Homebrew
brew install --cask miniconda
```

**Linux:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### HMMER Build Fails

```bash
# Check build logs
cd external/hmmer-3.4
./configure --prefix=$(pwd)/../hmmer
make

# If errors occur, install dependencies
# macOS:
brew install automake autoconf

# Linux:
sudo apt-get install build-essential
```

### CCMpred Build Fails (HPC)

Common issue: CUDA not found

```bash
# Check CUDA modules
module avail cuda

# Load appropriate version
module load cuda/11.0  # Or whatever is available

# Check CUDA is accessible
nvcc --version

# Rebuild
cd release
rm -rf external/CCMpred
bash setup_hpc.sh
```

### Disk Space Issues

The full installation requires ~5GB. Check available space:
```bash
df -h .
```

If space is limited, you can:
1. Use a different installation location with more space
2. Skip the UniRef90 database download (included in full repo but not this release)
3. Remove example data after testing: `rm -rf data/raw/*.raw` (saves 2GB)

### Python Package Installation Fails

If conda packages fail to install:

```bash
# Try with mamba (faster)
conda install -n base -c conda-forge mamba
mamba env create -f environment_mac.yml --prefix conda/

# Or specify conda-forge channel explicitly
conda env create -f environment_mac.yml --prefix conda/ -c conda-forge
```

## Clean Installation

To completely remove and reinstall:

```bash
# Remove conda environment and external tools
rm -rf conda/ external/

# Re-run setup
bash setup_mac.sh  # or setup_hpc.sh
```
