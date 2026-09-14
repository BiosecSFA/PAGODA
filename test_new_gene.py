#!/usr/bin/env python3
"""
Test New Gene with Conditionally Essential Genes

This script takes a new gene FASTA file, builds models for it, designs overlaps
with the 4 conditionally essential genes (purE, hisI, ilvE, gltA), and reports
the percentage of designs that meet Potts score thresholds for the experimentally
tested termini of the conditional essential genes.

Experimental configurations tested (only the validated gene orders):
  - purE C-terminus (upstream in purE_hisI library) + any downstream gene
  - hisI N-terminus (downstream in purE_hisI library) + any upstream gene
  - ilvE C-terminus (upstream in ilvE_gltA library) + any downstream gene
  - gltA N-terminus (downstream in ilvE_gltA library) + any upstream gene

Usage:
    python test_new_gene.py newgene.fasta --min-overlap 10 --max-overlap 50 --database data/db/uniref90.fasta
    python test_new_gene.py newgene.fasta --skip-model-gen  # Use if models already exist
"""

import argparse
import os
import sys
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

# Add repo to path
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

# Potts score thresholds determined from experimental enrichment data
# These are approximate values - run potts_score_threshold.py for exact values
POTTS_THRESHOLDS = {
    'purE': -108.10103607177734,   # From purE_hisI and hisI_purE analysis
    'hisI': -105.67639923095703,   # From purE_hisI and hisI_purE analysis
    'ilvE': -247.56973266601562,   # From ilvE_gltA analysis
    'gltA': -266.75701904296875,  # From ilvE_gltA analysis
}

# Experimentally tested termini for each cond ess gene
# Based on purE_hisI library: purE C-term (upstream) + hisI N-term (downstream)
TESTED_TERMINI = {
    'purE': {'C': True, 'N': False},   # Only C-terminus tested (purE upstream in purE_hisI)
    'hisI': {'C': False, 'N': True},   # Only N-terminus tested (hisI downstream in purE_hisI)
    'ilvE': {'C': True, 'N': False},   # Only C-terminus tested (ilvE upstream in ilvE_gltA)
    'gltA': {'C': False, 'N': True},   # Only N-terminus tested (gltA downstream in ilvE_gltA)
}


def is_experimentally_tested(cond_ess_gene, gene_position, order):
    """
    Check if this design tests an experimentally validated terminus.

    Args:
        cond_ess_gene: Name of conditional essential gene (purE/hisI/ilvE/gltA)
        gene_position: 1 or 2 (whether cond_ess gene is gene1 or gene2)
        order: "12" or "21"

    Returns:
        bool: True if this terminus was experimentally tested
    """
    # Determine which terminus of the cond ess gene is in the overlap
    if order == "12":
        # gene1 upstream (C-term), gene2 downstream (N-term)
        terminus = 'C' if gene_position == 1 else 'N'
    else:  # order == "21"
        # gene2 upstream (C-term), gene1 downstream (N-term)
        terminus = 'C' if gene_position == 2 else 'N'

    return TESTED_TERMINI[cond_ess_gene][terminus]


def parse_fasta(fasta_file):
    """Parse FASTA file and extract gene name and sequence."""
    record = list(SeqIO.parse(fasta_file, "fasta"))[0]
    gene_name = record.id
    protein_seq = str(record.seq)
    return gene_name, protein_seq



def prepare_gene_files(gene_name, protein_seq, fasta_file, data_dir="data", cds_file_path=None):
    """Prepare required input files for gene model generation."""
    # Create split_fasta entry
    split_fasta_file = f"{data_dir}/split_fasta/{gene_name}.fasta"
    os.makedirs(f"{data_dir}/split_fasta", exist_ok=True)

    with open(split_fasta_file, 'w') as f:
        f.write(f">{gene_name}\n{protein_seq}\n")
    print(f"  ✓ Created {split_fasta_file}")

    # Check if CDS file already exists
    split_cds_file = f"{data_dir}/split_cds/{gene_name}.fasta"
    os.makedirs(f"{data_dir}/split_cds", exist_ok=True)

    if os.path.exists(split_cds_file):
        print(f"  ✓ CDS file already exists: {split_cds_file}")
        return True

    # If CDS file provided via command line, copy it
    if cds_file_path:
        if not os.path.exists(cds_file_path):
            print(f"  ✗ CDS file not found: {cds_file_path}")
            return False

        # Read CDS from provided file
        cds_record = list(SeqIO.parse(cds_file_path, "fasta"))[0]
        cds_seq = str(cds_record.seq).upper()

        # Validate length
        expected_len = len(protein_seq) * 3
        if abs(len(cds_seq) - expected_len) > 10:
            print(f"  ⚠ Warning: CDS length ({len(cds_seq)} nt) doesn't match protein length ({len(protein_seq)} aa)")
            print(f"            Expected ~{expected_len} nt")

        # Write to split_cds
        with open(split_cds_file, 'w') as f:
            f.write(f">{gene_name}\n{cds_seq}\n")
        print(f"  ✓ Created {split_cds_file} from {cds_file_path}")
        return True

    # Prompt user for CDS sequence
    print(f"\n  CDS (nucleotide) sequence required for {gene_name}")
    print(f"  Please paste the coding sequence (DNA/RNA), then press Enter twice:")
    print(f"  (Paste multi-line sequences, then press Enter on empty line to finish)")
    print()

    cds_lines = []
    empty_count = 0

    while empty_count < 1:
        try:
            line = input()
            if line.strip() == "":
                empty_count += 1
            else:
                empty_count = 0
                cds_lines.append(line.strip())
        except EOFError:
            break

    if not cds_lines:
        print(f"\n  ✗ No CDS sequence provided")
        return False

    # Join and clean the sequence
    cds_seq = ''.join(cds_lines).upper()
    # Remove any whitespace, numbers, or common FASTA artifacts
    cds_seq = ''.join(c for c in cds_seq if c in 'ATGCUN')

    if len(cds_seq) == 0:
        print(f"\n  ✗ Invalid CDS sequence (no valid nucleotides)")
        return False

    # Validate: should be ~3x protein length
    expected_len = len(protein_seq) * 3
    if abs(len(cds_seq) - expected_len) > 10:
        print(f"\n  ⚠ Warning: CDS length ({len(cds_seq)} nt) doesn't match protein length ({len(protein_seq)} aa)")
        print(f"            Expected ~{expected_len} nt")
        response = input("  Continue anyway? [y/N] ").strip().lower()
        if response != 'y':
            return False

    # Write CDS file
    with open(split_cds_file, 'w') as f:
        f.write(f">{gene_name}\n{cds_seq}\n")
    print(f"  ✓ Created {split_cds_file}")

    return True


def run_model_generation(gene_name, database, data_dir="data", skip_ccmpred=False):
    """Run prepare_gene_models.py to generate MSA, HMM, and Potts models."""
    print(f"\n[Model Generation] {gene_name}")

    cmd = [
        "python", "scripts/prepare_gene_models.py",
        gene_name,
        "--database", database,
        "--data-dir", data_dir,
        "--skip-validation",
        "--force"
    ]

    if skip_ccmpred:
        cmd.append("--skip-ccmpred")
        print("  Note: Skipping CCMpred (Potts models will not be generated)")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ✗ Model generation failed:")
        if result.stdout:
            print("  STDOUT:")
            print(result.stdout)
        if result.stderr:
            print("  STDERR:")
            print(result.stderr)
        return False

    # Show successful output
    if result.stdout:
        print(result.stdout)

    # Check generated files (only HMM and Potts models are required)
    required_files = [
        f"{data_dir}/hmm/{gene_name}.hmm",
    ]

    if not skip_ccmpred:
        required_files.append(f"{data_dir}/raw/{gene_name}.raw")

    for f in required_files:
        if not os.path.exists(f):
            print(f"  ✗ Expected file not found: {f}")
            return False
        print(f"  ✓ {f}")

    return True


def run_dpalign(new_gene, cond_ess_gene, rep, min_overlap, max_overlap,
                output_dir, data_dir="data"):
    """Run dpalign.py to generate overlapping designs."""
    cmd = [
        "python", "scripts/dpalign.py",
        new_gene, cond_ess_gene, str(rep),
        "--del_penalty", "4.0",
        "--ins_penalty", "1.0",
        "--match_weight", "0.3",
        "--min_overlap", str(min_overlap),
        "--max_overlap", str(max_overlap),
        "--data_dir", data_dir,
        "--output_dir", output_dir
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"    ✗ dpalign failed:")
        print(result.stderr)
        return None

    # Find output file
    expected_file = f"{output_dir}/{new_gene}_{cond_ess_gene}_del4.0_ins1.0_0.3_rep{rep}.csv"

    if not os.path.exists(expected_file):
        print(f"    ✗ Output file not found: {expected_file}")
        return None

    return expected_file


def analyze_designs(csv_file, new_gene, cond_ess_gene):
    """
    Analyze designs and return percentage meeting Potts threshold
    for experimentally tested terminus of conditional essential gene.
    """
    df = pd.read_csv(csv_file)

    if len(df) == 0:
        return {
            'total_designs': 0,
            'tested_configs': 0,
            'passing_threshold': 0,
            'percentage': 0.0,
            'details': []
        }

    # Determine which gene is gene1 vs gene2
    is_new_gene_gene1 = (df['gene1'].iloc[0] == new_gene)
    cond_ess_position = 2 if is_new_gene_gene1 else 1
    cond_ess_score_col = f'optim{cond_ess_position}_score'

    threshold = POTTS_THRESHOLDS[cond_ess_gene]

    tested_configs = []
    for _, row in df.iterrows():
        order = str(row['order'])  # Convert to string for comparison
        is_tested = is_experimentally_tested(cond_ess_gene, cond_ess_position, order)

        if is_tested:
            score = row[cond_ess_score_col]
            passes = score >= threshold
            tested_configs.append({
                'order': order,
                'overlaplen': row['overlaplen'],
                'score': score,
                'passes': passes
            })

    n_tested = len(tested_configs)
    n_passing = sum(1 for c in tested_configs if c['passes'])
    percentage = (n_passing / n_tested * 100) if n_tested > 0 else 0.0

    return {
        'total_designs': len(df),
        'tested_configs': n_tested,
        'passing_threshold': n_passing,
        'percentage': percentage,
        'details': tested_configs
    }


def main():
    parser = argparse.ArgumentParser(
        description='Test new gene with conditionally essential genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Full model generation (MSA, HMM, Potts)
  python test_new_gene.py newgene.fasta --database data/db/uniref90.fasta

  # Generate MSA/HMM only (skip CCMpred)
  python test_new_gene.py newgene.fasta --database data/db/uniref90.fasta --skip-ccmpred

  # Non-interactive mode - provide CDS file
  python test_new_gene.py newgene.fasta --cds-file newgene_cds.fasta --skip-model-gen

  # Skip model generation (models already exist)
  python test_new_gene.py newgene.fasta --skip-model-gen

  # Custom overlap range
  python test_new_gene.py newgene.fasta --min-overlap 15 --max-overlap 40 --skip-model-gen
        '''
    )

    parser.add_argument('fasta', help='FASTA file with new gene protein sequence')
    parser.add_argument('--database', default='data/db/uniref90.fasta',
                       help='UniRef90 database for MSA generation')
    parser.add_argument('--min-overlap', type=int, default=10,
                       help='Minimum overlap length (aa, default: 10)')
    parser.add_argument('--max-overlap', type=int, default=30,
                       help='Maximum overlap length (aa, default: 50)')
    parser.add_argument('--data-dir', default='data',
                       help='Data directory (default: data)')
    parser.add_argument('--output-dir', default='test_new_gene_output',
                       help='Output directory (default: test_new_gene_output)')
    parser.add_argument('--skip-model-gen', action='store_true',
                       help='Skip model generation (models must already exist)')
    parser.add_argument('--skip-ccmpred', action='store_true', default=False,
                       help='Skip CCMpred step (default: False, generates Potts models)')
    parser.add_argument('--cds-file', type=str,
                       help='CDS FASTA file (optional, otherwise will prompt interactively)')

    args = parser.parse_args()

    # Parse input FASTA
    print("=" * 80)
    print("Testing New Gene with Conditionally Essential Genes")
    print("=" * 80)
    print(f"\n[Input] Parsing {args.fasta}")

    gene_name, protein_seq = parse_fasta(args.fasta)
    print(f"  Gene: {gene_name}")
    print(f"  Length: {len(protein_seq)} aa")
    print(f"  Sequence: {protein_seq[:50]}...")

    # Prepare gene files
    print(f"\n[Setup] Preparing gene files...")

    if not prepare_gene_files(gene_name, protein_seq, args.fasta, args.data_dir, args.cds_file):
        print("\n✗ Setup failed: CDS file required")
        print("  Exiting...")
        return 1

    # Generate models
    if not args.skip_model_gen:
        # Check for database file (uncompressed or compressed)
        if not os.path.exists(args.database):
            compressed_db = args.database + ".gz"
            if os.path.exists(compressed_db):
                print(f"\n[Database] Found compressed: {compressed_db}")
                print("  Decompressing (this may take a few minutes)...")
                import gzip
                import shutil
                with gzip.open(compressed_db, 'rb') as f_in:
                    with open(args.database, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                print(f"  ✓ Decompressed to {args.database}")
            else:
                print(f"\n✗ Database not found: {args.database}")
                print("  Download UniRef90 or use --skip-model-gen if models already exist")
                return 1

        success = run_model_generation(
            gene_name,
            args.database,
            args.data_dir,
            args.skip_ccmpred
        )

        if not success:
            print("\n✗ Model generation failed")
            return 1
    else:
        print(f"\n[Model Generation] Skipped (using existing models)")

    # Design overlaps with each conditional essential gene
    cond_ess_genes = ['purE', 'hisI', 'ilvE', 'gltA']

    print(f"\n[Design] Generating overlaps with conditional essential genes")
    print(f"  Overlap range: {args.min_overlap}-{args.max_overlap} aa")
    print(f"  Output directory: {args.output_dir}")

    os.makedirs(args.output_dir, exist_ok=True)

    results = {}

    for cond_ess_gene in cond_ess_genes:
        print(f"\n  {gene_name} × {cond_ess_gene}:")

        # Run dpalign
        csv_file = run_dpalign(
            gene_name, cond_ess_gene, 1,
            args.min_overlap, args.max_overlap,
            args.output_dir, args.data_dir
        )

        if csv_file is None:
            results[cond_ess_gene] = None
            continue

        # Analyze results
        analysis = analyze_designs(csv_file, gene_name, cond_ess_gene)
        results[cond_ess_gene] = analysis

        print(f"    Total designs: {analysis['total_designs']}")
        print(f"    Has reference gene: {analysis['tested_configs']}")
        print(f"    Passing Potts threshold (≥{POTTS_THRESHOLDS[cond_ess_gene]}): {analysis['passing_threshold']}")
        print(f"    Success rate: {analysis['percentage']:.1f}%")

    # Summary report
    print("\n" + "=" * 80)
    print("SUMMARY REPORT")
    print("=" * 80)
    print(f"\nNew gene: {gene_name} ({len(protein_seq)} aa)")
    print(f"Overlap range: {args.min_overlap}-{args.max_overlap} aa")
    print(f"\nPotts Score Thresholds:")
    for gene, thr in POTTS_THRESHOLDS.items():
        print(f"  {gene}: {thr}")

    print(f"\n{'Gene':<8} {'Tested':<10} {'Total':<10} {'Passing':<10} {'Success %':<12} {'Terminus'}")
    print("-" * 80)

    for cond_ess_gene in cond_ess_genes:
        if results[cond_ess_gene] is None:
            print(f"{cond_ess_gene:<8} {'FAILED':<10} {'-':<10} {'-':<10} {'-':<12} {'-'}")
            continue

        r = results[cond_ess_gene]
        tested_termini = [k for k, v in TESTED_TERMINI[cond_ess_gene].items() if v]
        terminus_str = '+'.join([f"{t}-term" for t in tested_termini])

        print(f"{cond_ess_gene:<8} "
              f"{r['tested_configs']:<10} "
              f"{r['total_designs']:<10} "
              f"{r['passing_threshold']:<10} "
              f"{r['percentage']:<12.1f} "
              f"{terminus_str}")

    print("\n" + "=" * 80)
    print("Test complete!")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
