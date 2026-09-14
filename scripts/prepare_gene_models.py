#!/usr/bin/env python3
"""
Prepare Gene Models Pipeline

This script takes a gene name and prepares all necessary files for overlapping gene design:
1. Extracts protein sequence from gene_info_final_table.tsv
2. Uses phmmer to search for homologs and build MSA
3. Trains HMM profile model using hmmbuild
4. Converts MSA to format for ccmpred
5. Trains Potts model using ccmpred

All output files are organized in the data/ directory structure.
"""

import argparse
import sys
import os
import subprocess
import tempfile
from pathlib import Path


def main():
    """Main function - imports are inside to speed up --help"""
    args = parse_arguments()
    
    # Import pandas only when needed
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    gene_name = args.gene
    data_dir = Path(args.data_dir)
    
    # File paths
    fasta_dir = data_dir / "split_fasta"
    cds_dir = data_dir / "split_cds"
    msa_dir = data_dir / "msa"
    hmm_dir = data_dir / "hmm"
    raw_dir = data_dir / "raw"
    
    # Create output directories
    for directory in [fasta_dir, cds_dir, msa_dir, hmm_dir, raw_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing gene: {gene_name}")
    print("=" * 60)

    # Step 1: Get protein sequence (try FASTA first, then table)
    print("\n[1/5] Loading protein sequence...")
    protein_seq = None
    uniprot_id = "unknown"

    # Try loading from existing FASTA file first
    fasta_file = fasta_dir / f"{gene_name}.fasta"
    if fasta_file.exists():
        print(f"  Found existing FASTA: {fasta_file}")
        record = list(SeqIO.parse(fasta_file, "fasta"))[0]
        protein_seq = str(record.seq)
        print(f"  Loaded protein sequence: {len(protein_seq)} amino acids")
    else:
        print(f"  ✗ Protein FASTA not found: {fasta_file}")
        print(f"  Please provide a protein sequence FASTA file at: {fasta_file}")
        print(f"  Format: >{gene_name}\\n<sequence>")
        sys.exit(1)

    # Step 1.5: Get DNA sequence (try split_cds first, then cds.fasta)
    print("\n  Loading DNA sequence...")

    cds_file = cds_dir / f"{gene_name}.fasta"
    cds_fasta_file = data_dir / "cds.fasta"
    dna_seq = None

    # Try split_cds first
    if cds_file.exists():
        print(f"  Found existing CDS: {cds_file}")
        record = list(SeqIO.parse(cds_file, "fasta"))[0]
        dna_seq = str(record.seq)
        print(f"  Loaded DNA sequence: {len(dna_seq)} nucleotides")

    # If not found, try cds.fasta
    elif cds_fasta_file.exists():
        print(f"  Searching in {cds_fasta_file}...")
        for record in SeqIO.parse(cds_fasta_file, "fasta"):
            if record.id == gene_name:
                dna_seq = str(record.seq)
                print(f"  ✓ Found DNA sequence: {len(dna_seq)} nucleotides")
                # Save to split_cds
                SeqIO.write(record, cds_file, "fasta")
                print(f"  Saved to: {cds_file}")
                break

        if dna_seq is None:
            print(f"  ✗ Gene {gene_name} not found in {cds_fasta_file}")

    if dna_seq is None:
        print(f"  ⚠ Warning: No DNA sequence found")
        print(f"  Provide CDS in: {cds_file}")
        print(f"  Or add to: {cds_fasta_file}")
        print(f"  DNA validation will be skipped")
    
    # Verify translation if we have DNA sequence
    if dna_seq:
        # Translate and verify
        translated = str(Seq(dna_seq).translate())
        
        # Remove stop codon if present
        if translated.endswith('*'):
            translated = translated[:-1]
        
        # Compare with protein sequence from table
        if translated == protein_seq:
            print(f"  ✓ DNA sequence translates correctly to protein sequence")
        else:
            print(f"  ✗ WARNING: DNA translation mismatch!")
            print(f"    Expected protein: {protein_seq[:50]}...")
            print(f"    Translated:       {translated[:50]}...")
            print(f"    Length - Expected: {len(protein_seq)}, Translated: {len(translated)}")
            
            # Try to align and show differences
            if len(translated) == len(protein_seq):
                diffs = sum(1 for a, b in zip(protein_seq, translated) if a != b)
                print(f"    Differences: {diffs} amino acids")
                
                # Show first few mismatches
                mismatches = [(i, a, b) for i, (a, b) in enumerate(zip(protein_seq, translated)) if a != b]
                if mismatches:
                    print(f"    First mismatches:")
                    for i, expected, got in mismatches[:5]:
                        print(f"      Position {i+1}: Expected {expected}, Got {got}")
            
            if not args.skip_validation:
                response = input("\n  Continue anyway? (y/n): ")
                if response.lower() != 'y':
                    sys.exit(1)
    else:
        print(f"  ✗ WARNING: No DNA sequence available for validation")
        print(f"  Continuing with protein sequence only...")

    # Protein FASTA already exists (we loaded it above), just reference it
    protein_fasta = fasta_file
    
    # Step 2: Run phmmer to generate MSA
    print("\n[2/5] Running phmmer to find homologs...")
    
    msa_file = msa_dir / f"{gene_name}.sto"
    
    # Check if MSA already exists
    if msa_file.exists() and not args.force:
        print(f"  ✓ MSA already exists: {msa_file}")
        # Count sequences in existing MSA
        with open(msa_file) as f:
            n_seqs = sum(1 for line in f if line.startswith('#=GS'))
        print(f"  Found {n_seqs} sequences in existing MSA")
        print("  Skipping phmmer search (use --force to regenerate)")
    else:
        # Check if phmmer is available
        try:
            subprocess.run(['phmmer', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Error: phmmer not found. Please install HMMER (≥3.4)")
            print("  macOS: brew install hmmer")
            print("  Linux: sudo apt-get install hmmer")
            sys.exit(1)
        
        # Use UniRef90 or specified database
        if args.database:
            db_path = Path(args.database)
            if not db_path.exists():
                print(f"Error: Database not found at {args.database}")
                sys.exit(1)
        else:
            print("  Note: No database specified. You need a protein sequence database.")
            print("  Use --database to specify path to UniRef90, UniProt, or similar FASTA database")
            print("  Skipping phmmer search...")
            
            # Create dummy MSA with just the query sequence for testing
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                tmp.write(f">{gene_name}\n{protein_seq}\n")
                tmp_fasta = tmp.name
            
            # Convert single sequence to Stockholm format
            subprocess.run([
                'esl-reformat', 'stockholm', tmp_fasta
            ], stdout=open(msa_file, 'w'), check=True)
            os.unlink(tmp_fasta)
            
            print(f"  Created single-sequence MSA: {msa_file}")
            db_path = None
        
        if db_path:
            # Calculate scaled bitscore threshold (EVcouplings default: 0.5 * sequence_length)
            seq_length = len(protein_seq)
            scaled_bitscore = args.bitscore_per_pos * seq_length
            
            # Run phmmer with EVcouplings-compatible settings
            phmmer_output = msa_dir / f"{gene_name}.phmmer"
            tblout = msa_dir / f"{gene_name}.tblout"
            
            cmd = [
                'phmmer',
                '--incT', str(scaled_bitscore),  # Inclusion threshold
                '-T', str(scaled_bitscore),  # Reporting threshold
                '--popen', '0.02',  # Gap open probability
                '--pextend', '0.4',  # Gap extend probability
                '--mx', 'BLOSUM62',  # Substitution matrix
                '--tblout', str(tblout),  # Table output
                '-A', str(msa_file),  # MSA output in Stockholm format
                '--noali',  # Don't include alignments in main output
                '--notextw',  # Unlimited text width
                '--cpu', str(args.cpu),  # Number of CPUs
                '-o', str(phmmer_output),  # Text output
                str(protein_fasta),
                str(db_path)
            ]
            
            print(f"  Scaled bitscore threshold: {scaled_bitscore:.1f} ({args.bitscore_per_pos} * {seq_length})")
            print(f"  Running: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"Error running phmmer: {result.stderr}")
                sys.exit(1)
            
            print(f"  MSA saved to: {msa_file}")
            print(f"  Phmmer output: {phmmer_output}")
            
            # Count sequences in MSA
            with open(msa_file) as f:
                n_seqs = sum(1 for line in f if line.startswith('#=GS'))
            print(f"  Found {n_seqs} homologous sequences")
    
    # Step 3: Build HMM profile
    print("\n[3/5] Building HMM profile with hmmbuild...")
    
    hmm_file = hmm_dir / f"{gene_name}.hmm"
    
    # Check if HMM already exists
    if hmm_file.exists() and not args.force:
        print(f"  ✓ HMM profile already exists: {hmm_file}")
        print("  Skipping hmmbuild (use --force to regenerate)")
    else:
        cmd = [
            'hmmbuild',
            '--amino',  # Amino acid sequences
            '--hand',   # Manual alignment, do not realign
            '--cpu', str(args.cpu),  # Number of CPUs
            str(hmm_file),
            str(msa_file)
        ]
        
        print(f"  Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running hmmbuild: {result.stderr}")
            sys.exit(1)
        
        print(f"  HMM profile saved to: {hmm_file}")
        
        # Press HMM for faster searches (creates .h3m, .h3p, .h3i, .h3f files)
        print("  Pressing HMM for indexing...")
        cmd = ['hmmpress', '-f', str(hmm_file)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Warning: hmmpress failed: {result.stderr}")
        else:
            print("  ✓ HMM indexed successfully")
    
    # Step 4: Convert MSA to aligned format for ccmpred
    print("\n[4/5] Converting MSA to aligned format for ccmpred...")
    
    aln_file = msa_dir / f"{gene_name}.aln"
    
    # Check if alignment file already exists
    if aln_file.exists() and not args.force:
        print(f"  ✓ Alignment file already exists: {aln_file}")
        print("  Skipping conversion (use --force to regenerate)")
    else:
        # Convert Stockholm to .aln format, keeping only reference sequence positions
        try:
            from Bio import AlignIO
            
            # Read the Stockholm alignment
            msa_alignment = AlignIO.read(msa_file, "stockholm")
            
            # Get reference annotation (marks positions in the query/reference sequence)
            if 'reference_annotation' in msa_alignment.column_annotations:
                ref_annotation = msa_alignment.column_annotations['reference_annotation']
            else:
                # If no reference annotation, use first sequence positions
                ref_seq = str(msa_alignment[0].seq)
                ref_annotation = ''.join('x' if c != '-' and c != '.' else '.' for c in ref_seq)
            
            # Identify positions where reference has a residue (not a gap)
            ref_positions = [i for i, char in enumerate(ref_annotation) if char != '.']
            
            print(f"  Reference sequence length: {len(ref_positions)} positions (out of {msa_alignment.get_alignment_length()} total)")
            
            # Write .aln file with only reference positions (no gaps in reference)
            with open(aln_file, 'w') as f:
                for record in msa_alignment:
                    seq = str(record.seq)
                    # Extract only positions where reference has residues
                    filtered_seq = ''.join(seq[i] for i in ref_positions)
                    f.write(filtered_seq + '\n')
            
            print(f"  Alignment file saved to: {aln_file}")
            print(f"  Filtered to {len(ref_positions)} positions matching reference sequence")
            
        except Exception as e:
            print(f"Error converting to .aln format: {e}")
            import traceback
            traceback.print_exc()
            print("  Continuing without alignment conversion...")
            aln_file = None
    
    # Step 5: Run ccmpred to train Potts model
    if args.skip_ccmpred:
        print("\n[5/5] Skipping ccmpred (--skip-ccmpred flag set)")
        print("  MSA alignment ready for ccmpred: {}".format(aln_file if aln_file else "N/A"))
        print("  Transfer .aln file to GPU cluster and run ccmpred separately")
        sys.exit()
    else:
        print("\n[5/5] Training Potts model with ccmpred...")
        
        # Check if ccmpred is available (try PATH first, then external/CCMpred/bin)
        ccmpred_cmd = 'ccmpred'
        ccmpred_available = False
    
    # Try ccmpred in PATH
    try:
        result = subprocess.run(['ccmpred', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # ccmpred -h may return non-zero, so just check if it runs without FileNotFoundError
        ccmpred_available = True
    except FileNotFoundError:
        pass
    
    # If not in PATH, try external/CCMpred/bin/ccmpred relative to current working directory
    if not ccmpred_available:
        # Try relative to current working directory first
        ccmpred_local = Path.cwd() / 'external' / 'CCMpred' / 'bin' / 'ccmpred'
        print(ccmpred_local)
        if not ccmpred_local.exists():
            # Fall back to relative to script location
            script_dir = Path(__file__).parent.parent
            ccmpred_local = script_dir / 'external' / 'CCMpred' / 'bin' / 'ccmpred'
        
        if ccmpred_local.exists():
            try:
                result = subprocess.run([str(ccmpred_local), '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                ccmpred_cmd = str(ccmpred_local)
                ccmpred_available = True
                print(f"  Using ccmpred from: {ccmpred_local}")
            except (FileNotFoundError, PermissionError) as e:
                print(f"  Warning: Found ccmpred at {ccmpred_local} but cannot execute: {e}")
    
    if not ccmpred_available:
        print("")
        print("=" * 70)
        print("WARNING: ccmpred not found in PATH")
        print("=" * 70)
        print("")
        print("  CCMpred is required to train Potts models (.raw files).")
        print("  The pipeline REQUIRES these Potts models to function.")
        print("")
        print("  Options:")
        print("  1. If setup_environment.sh installed a CPU-only version:")
        print("     - Check: external/ccmpred/bin/ccmpred")
        print("     - Run: source activate_pipeline.sh")
        print("")
        print("  2. Install CCMpred on a GPU system for better performance:")
        print("     - Transfer MSA files (data/msa/*.aln) to GPU system")
        print("     - Install: https://github.com/soedinglab/CCMpred")
        print("     - Run: ccmpred {}.aln {}.raw".format(gene_name, gene_name))
        print("     - Transfer .raw files back to data/raw/")
        print("")
        print("  3. Use pre-existing .raw files if available")
        print("")
        print("=" * 70)
        
        # Check if we can continue anyway
        raw_file = raw_dir / f"{gene_name}.raw"
        if raw_file.exists():
            print(f"\n✓ Found existing Potts model: {raw_file}")
            print("  Skipping ccmpred training step")
        else:
            print("\nPipeline cannot continue without Potts models.")
            print("Please install ccmpred or obtain .raw files from a GPU system.")
            sys.exit(1)
    else:
        # CCMpred is available, run it
        if aln_file:
            mat_file = raw_dir / f"{gene_name}.mat"
            raw_file = raw_dir / f"{gene_name}.raw"
            
            cmd = [
                ccmpred_cmd,
                "-r",  str(raw_file),
                str(aln_file),
                str(mat_file),
            ]
            
            # Add iteration count (default: 100)
            cmd.extend(['-n', str(args.ccmpred_niter)])
            
            # Add thread count if specified
            if args.ccmpred_threads:
                cmd.extend(['-t', str(args.ccmpred_threads)])
            
            # Output log file
            ccmpred_log = raw_dir / f"{gene_name}.ccmpred.log"
            
            print(f"  Running: {' '.join(cmd)}")
            print(f"  CCMpred iterations: {args.ccmpred_niter}")
            print(f"  Progress log: {ccmpred_log}")
            print("  Note: This may take several minutes depending on sequence length...")
            
            # Run ccmpred and redirect stderr to log file (unbuffered)
            with open(ccmpred_log, 'w', buffering=1) as log_file:
                result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, text=True, bufsize=1)
            
            if result.returncode != 0:
                print(f"Error running ccmpred. Check log: {ccmpred_log}")
                # Print last few lines of log
                with open(ccmpred_log, 'r') as f:
                    lines = f.readlines()
                    print("Last 10 lines of log:")
                    print(''.join(lines[-10:]))
                sys.exit(1)
            else:
                print(f"  ✓ Potts model saved to: {raw_file}")
                print(f"  ✓ Log saved to: {ccmpred_log}")
    
    # Summary
    print("\n" + "=" * 60)
    print("Pipeline complete! Generated files:")
    print(f"  Protein FASTA: {protein_fasta}")
    if msa_file.exists():
        print(f"  MSA (Stockholm): {msa_file}")
    if aln_file and aln_file.exists():
        print(f"  MSA (aligned): {aln_file}")
    if hmm_file.exists():
        print(f"  HMM profile: {hmm_file}")
    raw_file = raw_dir / f"{gene_name}.raw"
    if raw_file.exists():
        print(f"  Potts model: {raw_file}")
    
    print("\nReady for overlapping gene design!")


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Prepare all necessary models for a gene (MSA, HMM, Potts)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s trpC
  %(prog)s argH --database /path/to/uniref90.fasta
  %(prog)s hisF --ccmpred-niter 150 --ccmpred-threads 8
  %(prog)s purE --bitscore-per-pos 0.3 --force

Requirements:
  - HMMER (≥3.4): phmmer, hmmbuild, hmmpress, esl-reformat
  - ccmpred: for Potts model training
  - Protein sequence database (e.g., UniRef90) for homolog search

MSA Generation (EVcouplings-compatible):
  Uses phmmer with scaled bitscore thresholds.
  Default settings: 0.5 bits/position, BLOSUM62 matrix.
  Gap penalties: --popen 0.02, --pextend 0.4

Output structure:
  data/split_fasta/{gene}.fasta  - Query protein sequence
  data/msa/{gene}.sto            - Multiple sequence alignment (Stockholm)
  data/msa/{gene}.aln            - MSA in aligned format for ccmpred
  data/hmm/{gene}.hmm            - Profile HMM
  data/raw/{gene}.raw            - Potts model (ccmpred output)
        '''
    )
    
    parser.add_argument('gene', type=str,
                       help='Gene name (must exist in gene_info_final_table.tsv)')
    parser.add_argument('--database', type=str, default='data/db/uniref90.fasta',
                       help='Path to protein sequence database for phmmer search (default: data/db/uniref90.fasta)')
    parser.add_argument('--data-dir', type=str, default='data',
                       help='Base data directory (default: data)')
    parser.add_argument('--bitscore-per-pos', type=float, default=0.5,
                       help='Bitscore threshold per position for phmmer (default: 0.5, EVcouplings default)')
    parser.add_argument('--cpu', type=int, default=1,
                       help='Number of CPUs for phmmer (default: 1)')
    parser.add_argument('--ccmpred-threads', type=int,
                       help='Number of threads for ccmpred (default: auto)')
    parser.add_argument('--ccmpred-niter', type=int, default=250,
                       help='Number of ccmpred iterations (default: 250'
                            ')')
    parser.add_argument('--skip-validation', action='store_true',
                       help='Skip DNA translation validation prompts')
    parser.add_argument('--skip-ccmpred', action='store_true',
                       help='Skip ccmpred step (run MSA/HMM only, useful for running ccmpred separately on GPU cluster)')
    parser.add_argument('--force', action='store_true',
                       help='Force regeneration of existing files (MSA, HMM, etc.)')
    
    return parser.parse_args()


if __name__ == "__main__":
    main()
