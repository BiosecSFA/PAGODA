#!/bin/bash
# Activation script for Overlapping Gene Design Pipeline Release

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Initialize conda for this shell session
eval "$(conda shell.bash hook 2>/dev/null)" || eval "$(conda shell.zsh hook 2>/dev/null)" || true

# Activate conda environment from local conda/ directory
source activate "$SCRIPT_DIR/conda"

# CRITICAL: Disable user site-packages to prevent conflicts
# This prevents Python from loading packages from ~/.local/lib/python*/site-packages
export PYTHONNOUSERSITE=1

# Add HMMER and Easel to PATH (from conda/external/)
if [ -d "$SCRIPT_DIR/conda/external/hmmer/bin" ]; then
    export PATH="$SCRIPT_DIR/conda/external/hmmer/bin:$PATH"
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
echo "CCMpred: Not available on macOS (use pre-computed .raw files)"
echo ""
echo "Example commands:"
echo ""
echo "  # Design overlaps with controlled overlap range (10-30 aa)"
echo "  python scripts/dpalign.py --min_overlap 10 --max_overlap 30 purE hisI 1"
echo ""
python test_new_gene.py data/split_fasta/infA.fasta --database /p/lustre1/xu26/datangle/dbs/uniref90.fasta.gz --cds-file data/split_cds/infA.fasta 

echo "  # Test new gene against conditionally essential genes"
echo "  python test_new_gene.py data/split_fasta/infA.fasta \\"
echo "      --database /p/lustre1/xu26/datangle/dbs/uniref90.fasta.gz \\"
echo "      --cds-file data/split_cds/infA.fasta \\"
echo "      --min-overlap 10 --max-overlap 30"
echo ""
echo "  # Visualize Potts model couplings"
echo "  python scripts/plot_ccmpred_circos.py data/raw/purE.raw purE_circos.png \\"
echo "      --fasta data/split_fasta/purE.fasta"
echo ""
echo "To run tests:"
echo "  bash test_pipeline.sh"
echo ""
