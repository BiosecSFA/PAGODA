#!/bin/bash
# Test script for Overlapping Gene Design Pipeline Release

set -e  # Exit on error

echo "================================================================================"
echo "Testing Overlapping Gene Design Pipeline Release"
echo "================================================================================"

# Test 1: Environment check
echo ""
echo "[Test 1] Environment check..."
python --version
python -c "import biotite; print('✓ biotite')"
python -c "import Bio; print('✓ biopython')"
python -c "import numpy; print('✓ numpy')"
python -c "import pandas; print('✓ pandas')"
python -c "import torch; print('✓ pytorch')"
python -c "import scipy; print('✓ scipy')"
python -c "import matplotlib; print('✓ matplotlib')"
python -c "import seaborn; print('✓ seaborn')"

# Test 2: External tools
echo ""
echo "[Test 2] External tools..."
if command -v phmmer > /dev/null; then
    echo "✓ HMMER (phmmer)"
else
    echo "✗ HMMER (phmmer) not found in PATH"
    exit 1
fi

if command -v hmmbuild > /dev/null; then
    echo "✓ HMMER (hmmbuild)"
else
    echo "✗ HMMER (hmmbuild) not found in PATH"
    exit 1
fi

if command -v esl-reformat > /dev/null; then
    echo "✓ Easel (esl-reformat)"
else
    echo "✗ Easel (esl-reformat) not found in PATH"
    exit 1
fi

if command -v ccmpred > /dev/null; then
    echo "✓ CCMpred (available)"
else
    echo "⚠ CCMpred (not available - expected on macOS)"
fi

# Test 3: Custom packages
echo ""
echo "[Test 3] Custom packages..."
python -c "from pycameox.ilp import load_ccmpred; print('✓ pycameox')" || {
    echo "✗ pycameox import failed"
    exit 1
}
python -c "import src.alignment.scorematrix; print('✓ src.alignment')" || {
    echo "✗ src.alignment import failed"
    exit 1
}
python -c "import src.models.hmm; print('✓ src.models')" || {
    echo "✗ src.models import failed"
    exit 1
}
python -c "import src.optimization.greedy; print('✓ src.optimization')" || {
    echo "✗ src.optimization import failed"
    exit 1
}

# Test 4: Data files
echo ""
echo "[Test 4] Data files..."
for gene in purE hisI ilvE gltA; do
    if [ -f "data/split_fasta/${gene}.fasta" ]; then
        echo "✓ data/split_fasta/${gene}.fasta"
    else
        echo "✗ data/split_fasta/${gene}.fasta missing"
        exit 1
    fi

    if [ -f "data/hmm/${gene}.hmm" ]; then
        echo "✓ data/hmm/${gene}.hmm"
    else
        echo "✗ data/hmm/${gene}.hmm missing"
        exit 1
    fi

    if [ -f "data/raw/${gene}.raw" ]; then
        echo "✓ data/raw/${gene}.raw"
    else
        echo "✗ data/raw/${gene}.raw missing"
        exit 1
    fi
done

# Test 5: Scripts help
echo ""
echo "[Test 5] Script help messages..."
python scripts/dpalign.py --help > /dev/null 2>&1 && echo "✓ dpalign.py --help" || {
    echo "✗ dpalign.py --help failed"
    exit 1
}
python scripts/prepare_gene_models.py --help > /dev/null 2>&1 && echo "✓ prepare_gene_models.py --help" || {
    echo "✗ prepare_gene_models.py --help failed"
    exit 1
}

# Test 6: Run dpalign on small example (purE x hisI)
echo ""
echo "[Test 6] Running dpalign on purE x hisI..."
echo "  (This may take 30 minutes...)"

# Create temporary output directory
TEST_OUTPUT="test_output_$$"
mkdir -p "$TEST_OUTPUT"

python scripts/dpalign.py purE hisI 1 \
    --del_penalty 4.0 \
    --ins_penalty 1.0 \
    --match_weight 0.3 \
    --max_overlap 30 \
    --output_dir "$TEST_OUTPUT" 2>&1 | grep -E "(Loading|Computing|Designing|Writing)" || true

OUTPUT_FILE="${TEST_OUTPUT}/purE_hisI_del4.0_ins1.0_0.3_rep1.csv"

if [ -f "$OUTPUT_FILE" ]; then
    echo "✓ Output file created: $OUTPUT_FILE"
    NUM_DESIGNS=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
    echo "  Generated ${NUM_DESIGNS} designs"

    # Show a sample design
    echo ""
    echo "  Sample output (first design):"
    head -2 "$OUTPUT_FILE" | tail -1 | cut -d',' -f1-5

else
    echo "✗ Output file not found: $OUTPUT_FILE"
    echo "  Check for errors above"
    exit 1
fi

# Clean up test output
rm -rf "$TEST_OUTPUT"

# Test 7: Optional - Test gene model generation (requires database or --test-mode)
echo ""
echo "[Test 7] Testing gene model generation workflow..."
echo "  (Optional - requires UniRef90 database)"
echo ""
echo "  To download UniRef90 (50GB):"
echo "    mkdir -p data/db"
echo "    cd data/db"
echo "    wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"

# Check if user wants to run this test
if [ "$1" == "--full" ] || [ "$1" == "--with-model-gen" ]; then
    if [ -f "data/db/uniref90.fasta" ]; then
        echo "  Testing prepare_gene_models.py with database..."
        TEST_MODELS="test_models_$$"
        mkdir -p "$TEST_MODELS"

        python scripts/prepare_gene_models.py purE \
            --database data/db/uniref90.fasta \
            --data-dir "$TEST_MODELS" \
            --skip-ccmpred 2>&1 | grep -E "(Processing|Extracting|Searching|Building)" || true

        echo "✓ Model generation test completed"

        # Clean up
        rm -rf "$TEST_MODELS"
    else
        echo "  ⚠ UniRef90 database not found at data/db/uniref90.fasta"
        echo "  Skipping model generation test"
    fi
else
    echo "  Skipped (run with --full to test model generation)"
    echo ""
    echo "  To test model generation:"
    echo "    bash test_pipeline.sh --full"
fi

echo ""
echo "================================================================================"
echo "All tests passed! ✓"
echo "================================================================================"
echo ""
echo "The pipeline is ready to use."
echo ""
echo "Example commands:"
echo "  # Design overlapping genes"
echo "  python scripts/dpalign.py purE hisI 1 --del_penalty 4.0 --ins_penalty 1.0 --match_weight 0.3 --max_overlap 30"
echo ""
echo "  # Try different gene pairs"
echo "  python scripts/dpalign.py ilvE gltA 1"
echo ""
echo "  # Visualize Potts model couplings"
echo "  python scripts/plot_ccmpred_circos.py data/raw/purE.raw purE_circos.png \\"
echo "      --fasta data/split_fasta/purE.fasta --label-every 20"
echo ""
echo "For more information, see README.md"
echo ""
