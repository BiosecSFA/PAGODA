# Overlapping Gene Design Pipeline - Release Package

A computational pipeline for designing synthetic overlapping gene pairs using dynamic programming, Hidden Markov Models (HMMs), and Potts statistical models.

## Overview

This pipeline designs overlapping gene sequences where two coding sequences share a common nucleotide region while maintaining proper translation in different reading frames. The approach combines:

- **Dynamic Programming**: Optimal alignment with customizable gap penalties
- **Profile HMMs**: Evolutionary constraints from multiple sequence alignments
- **Potts Models**: Direct coupling analysis for sequence optimization
- **Greedy Optimization**: Codon-level refinement for improved fitness

## Quick Start

### macOS
```bash
bash setup_mac.sh
source activate_release.sh
bash test_pipeline.sh
python scripts/dpalign.py purE hisI 1
```

### HPC/Linux
```bash
module load cuda  # Load CUDA module first
bash setup_hpc.sh
source activate_release.sh
bash test_pipeline.sh
python scripts/dpalign.py purE hisI 1
```

## Requirements

- **Operating System**: macOS or Linux (HPC)
- **Package Manager**: Conda or Mamba (system-level installation required for initial setup)
- **Python**: 3.11 (installed locally in conda/)
- **CUDA** (HPC only): For Potts model training with CCMpred
- **Disk Space**: ~5GB for complete installation
  - 2GB: Example gene data (Potts models)
  - 3GB: Conda environment with all tools

**Fully Self-Contained**: All packages, tools, and dependencies are installed within the `conda/` directory. No global environment pollution. The entire environment can be transferred by copying just the `conda/` directory.

## Installation
### cloning the repo
```bash
git lfs install
git clone https://github.com/BiosecSFA/PAGODA.git
```
### macOS Setup

```bash
# 1. Install Miniconda if not already installed
# Download from: https://docs.conda.io/en/latest/miniconda.html

# 2. Run setup script
bash setup_mac.sh

# 3. Activate environment
source activate_release.sh
```

**Note**: CCMpred on macOS runs in CPU-only mode (slower than GPU). For faster Potts model generation, use HPC with CUDA. Pre-computed Potts models are included for the example genes.

### HPC/Linux Setup

```bash
# 1. Load CUDA module (if available)
module load cuda

# 2. Run setup script
bash setup_hpc.sh

# 3. Activate environment
source activate_release.sh
```

The setup script will:
- Create self-contained conda environment in `conda/` directory
- Install Python 3.11 and all dependencies locally
- Build HMMER 3.4 from source in `conda/external/hmmer/`
- Build CCMpred in `conda/external/ccmpred/` (CPU-only on macOS, CUDA on HPC)
- Install pycameox and potts packages in `conda/external/`
- Configure environment activation script

**Note**: Everything is installed within the `conda/` directory, making it fully portable and self-contained. Simply move or copy the `conda/` folder to transfer the entire environment.

## Testing

Run the test script to validate the installation:

```bash
bash test_pipeline.sh
```

This will verify:
- All Python packages are installed correctly
- External tools (HMMER, CCMpred) are available
- Data files exist for example genes
- Scripts can be executed
- Pipeline produces output for a test case

### Full Test (including model generation)

```bash
# Requires UniRef90 database (see "Generating Models" below)
bash test_pipeline.sh --full
```

## Example Genes

This release includes 4 example genes with complete data files:

| Gene | Size | Description |
|------|------|-------------|
| **purE** | 169 aa | Small, used in experimental validation (purE_hisI library) |
| **hisI** | 202 aa | Small-medium, used in experimental validation |
| **ilvE** | 309 aa | Medium, used in experimental validation (ilvE_gltA library) |
| **gltA** | 427 aa | Medium-large, used in experimental validation |

Each gene includes:
- Protein sequence (FASTA): `data/split_fasta/{gene}.fasta`
- CDS sequence (nucleotide): `data/split_cds/{gene}.fasta` OR `data/cds.fasta`
- Profile HMM: `data/hmm/{gene}.hmm`
- Potts model: `data/raw/{gene}.raw`

**Note on CDS sequences**: The pipeline requires individual CDS files in `data/split_cds/`. You can provide these in two ways:
1. **Individual files**: Place each gene's CDS sequence in `data/split_cds/{gene}.fasta`
2. **Combined file**: Place all CDS sequences in `data/cds.fasta` and extract them:
   ```bash
   # Extract individual CDS files from combined FASTA
   cd data
   for gene in gene1 gene2 gene3; do
     grep -A1 "^>${gene}$" cds.fasta > split_cds/${gene}.fasta
   done
   ```

## Running Examples

### Design Overlapping Genes

```bash
# Basic usage
python scripts/dpalign.py purE hisI 1

# With custom parameters
python scripts/dpalign.py ilvE gltA 1 \
    --del_penalty 4.0 \
    --ins_penalty 1.0 \
    --match_weight 0.3 \
    --output-dir results
```

**Parameters:**
- `gene1 gene2`: Names of genes to overlap
- `rep_number`: Replicate number (for tracking multiple runs)
- `--del_penalty`: Deletion penalty (higher = fewer deletions, default: 4.0)
- `--ins_penalty`: Insertion penalty (higher = fewer insertions, default: 1.0)
- `--match_weight`: Wild-type sequence preference (0-1, default: 0.3)
- `--output-dir`: Output directory (default: current directory)

**Output:**
- CSV file: `{output_dir}/{gene1}_{gene2}_del{pdel}_ins{pins}_{wmatch}_rep{rep}.csv`
- Contains all designed overlapping sequences with scores and properties

### Visualize Potts Model

```bash
# Generate circos plot showing residue couplings
python scripts/plot_ccmpred_circos.py \
    data/raw/purE.raw \
    purE_circos.png \
    --fasta data/split_fasta/purE.fasta \
    --label-every 20
```

## Generating Models for New Genes

To prepare models for genes not included in this release, you'll need the UniRef90 database (~50GB).

### Download UniRef90 Database

```bash
mkdir -p data/db
cd data/db
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
cd ../..
```

### Prepare Your Gene Data

Before generating models, ensure you have the required input files:

1. **Add gene information** to `data/gene_info_final_table.tsv`:
   ```
   Gene    Description    Length_aa
   newgene Description    300
   ```

2. **Add protein sequence** to `data/split_fasta/newgene.fasta`:
   ```
   >newgene
   MKKL...
   ```

3. **Add CDS sequence** - Choose one method:
   
   **Method A**: Individual file in `data/split_cds/newgene.fasta`:
   ```
   >newgene
   ATGAAAAAA...
   ```
   
   **Method B**: Add to combined file `data/cds.fasta`, then extract:
   ```bash
   cd data
   grep -A1 "^>newgene$" cds.fasta > split_cds/newgene.fasta
   ```

### Generate Models

**On macOS** (complete models, CCMpred in CPU-only mode):
```bash
python scripts/prepare_gene_models.py newgene \
    --database data/db/uniref90.fasta
```

Note: CCMpred runs slower on macOS CPU. For large proteins or many genes, use HPC with CUDA.

**On HPC/Linux** (complete models):
```bash
module load cuda
python scripts/prepare_gene_models.py newgene \
    --database data/db/uniref90.fasta
```

This will create:
1. MSA (multiple sequence alignment) in data/msa/
2. Profile HMM in data/hmm/
3. Potts model in data/raw/ (HPC only)

## Output Format

The pipeline generates CSV files with one row per designed overlapping sequence:

| Column | Description |
|--------|-------------|
| `seqid` | Unique sequence identifier |
| `gene1`, `gene2` | Gene names |
| `minlen` | Minimum overlap length constraint |
| `order` | Gene order (12 or 21) |
| `full_seq` | Complete overlapping nucleotide sequence |
| `start1`, `end1` | Coordinates in gene1 protein (nucleotide positions) |
| `start2`, `end2` | Coordinates in gene2 protein (nucleotide positions) |
| `overlaplen` | Actual overlap length (amino acids) |
| `raw1`, `raw2` | Initial aligned protein sequences (with gaps `-`) |
| `init1`, `init2` | Initial protein sequences (no gaps) |
| `optim1`, `optim2` | Optimized protein sequences (after greedy refinement) |
| `init1_score`, `init2_score` | Initial Potts model scores |
| `optim1_score`, `optim2_score` | Optimized Potts model scores |

## Troubleshooting

### macOS: CCMpred Not Available

CCMpred requires CUDA which is only available for NVIDIA GPUs. macOS does not support NVIDIA GPUs.

**Solutions:**
1. Use pre-computed Potts models (included for example genes)
2. Generate MSA and HMM on macOS, then run CCMpred on HPC:
   ```bash
   # On macOS
   python scripts/prepare_gene_models.py gene --database data/db/uniref90.fasta --skip-ccmpred
   
   # Transfer data/msa/gene.a3m to HPC
   
   # On HPC with GPU
   module load cuda
   ccmpred data/msa/gene.a3m data/raw/gene.raw
   
   # Transfer data/raw/gene.raw back to macOS
   ```

### HPC: CUDA Not Found

If CCMpred compilation fails:

```bash
# Load CUDA module
module avail cuda  # Check available versions
module load cuda/11.0  # Or appropriate version

# Rebuild CCMpred
rm -rf external/CCMpred
bash setup_hpc.sh
```

### Test Failures

Run the test script with verbose output:
```bash
bash -x test_pipeline.sh
```

Common issues:
- **Import errors**: Check that PYTHONPATH is set correctly in activate_release.sh
- **HMMER not found**: Verify external/hmmer/bin/ is in PATH
- **Missing data files**: Re-run setup script to rebuild

### Environment Activation Issues

If `source activate_release.sh` fails:
```bash
# Manually activate the local conda environment
conda activate $(pwd)/conda

# Set paths manually
export PATH="$(pwd)/external/hmmer/bin:$PATH"
export PYTHONPATH="$(pwd):$(pwd)/external/pycameox/src:$(pwd)/external/potts:$PYTHONPATH"
```

### Moving the Release Directory

The environment is fully self-contained and can be moved to a different location:
```bash
# Move the entire directory
mv release /new/location/

# Re-activate (paths are relative)
cd /new/location/release
source activate_release.sh
```

All paths in `activate_release.sh` are relative to the release directory, so it remains portable.

### Transferring Just the Environment

The `conda/` directory contains everything needed:
```bash
# Transfer to another system
tar -czf overlap-design-env.tar.gz conda/
# Transfer file, then extract on target system
tar -xzf overlap-design-env.tar.gz

# The conda/ directory now contains:
# - Python 3.11 and all packages
# - HMMER, pycameox, potts
# - CCMpred (if built on HPC)
```

## Algorithm Overview

### Dynamic Programming Matrix

The score matrix `S[i,j,k,l]` represents:
- `i`: Position in protein 1
- `j`: Position in protein 2
- `k`: Previous nucleotide (0-4, or deletion state)
- `l`: Current nucleotide (0-3)

States: `m` (match), `d` (deletion), `i1` (insertion seq1), `i2` (insertion seq2)

### Scoring Function

Total score = HMM log-probability + Potts energy

**HMM component:**
- Match emission probabilities from profile HMM
- Transition probabilities with adjustable penalties
- Insertion/deletion costs

**Potts component:**
- Single-site terms h(i, a)
- Pairwise coupling terms W(i, a, j, b)
- Wild-type sequence bonus (controlled by `match_weight`)

### Optimization Strategy

1. **Initial Design**: DP algorithm with HMM constraints finds optimal alignment
2. **Greedy Refinement**: For each codon in overlap, try all synonymous codons and keep changes that improve Potts score while respecting both reading frames

## Citation

If you use this pipeline, please cite:

[Publication information to be added]

## References

- **HMMER**: http://hmmer.org/
- **CCMpred**: https://github.com/soedinglab/CCMpred
- **Potts Models**: Weigt et al., PNAS 2009; Morcos et al., PNAS 2011
- **Profile HMMs**: Eddy, Current Opinion in Structural Biology 1996
- **pyCameoX**: Codon optimization library (BiosecSFA)
- **Biotite**: Kunzmann & Hamacher, BMC Bioinformatics 2018

## License

MIT License

## Contact
Chenling Xu: chenlingantelope@gmail.com

## Acknowledgments

This release package includes example genes with experimental validation data from the ilvE/gltA and purE_hisI libraries.
