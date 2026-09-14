#!/usr/bin/env python3
"""
Regenerate MSA for multiple genes with higher bitscore threshold in parallel

This script submits sbatch jobs to regenerate MSAs with a higher bitscore-per-pos
to generate smaller MSAs that are more manageable for ccmpred on GPU.
"""

import subprocess
import sys
import time
import re
from pathlib import Path


def submit_msa_job(gene, bitscore_per_pos=0.7):
    """Submit sbatch job to regenerate MSA for a single gene"""
    cmd = [
        'sbatch',
        f'--job-name=msa_{gene}',
        'batch/regenerate_msa.sh',
        gene,
        str(bitscore_per_pos)
    ]

    print(f"Submitting job for {gene}...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: Failed to submit job for {gene}")
        print(result.stderr)
        return None

    # Extract job ID from output (format: "Submitted batch job 12345")
    match = re.search(r'Submitted batch job (\d+)', result.stdout)
    if match:
        job_id = match.group(1)
        print(f"✓ Submitted job {job_id} for {gene}")
        return job_id
    else:
        print(f"WARNING: Could not parse job ID from: {result.stdout}")
        return None


def check_job_status(job_ids):
    """Check status of submitted jobs"""
    if not job_ids:
        return {}

    cmd = ['squeue', '--job', ','.join(job_ids), '--format=%i %T', '--noheader']
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        return {}

    status_map = {}
    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split()
            if len(parts) >= 2:
                status_map[parts[0]] = parts[1]

    return status_map


def main():
    """Submit MSA regeneration jobs for multiple genes in parallel"""
    genes = ['metR', 'cysB', 'trpD', 'trpE', 'cysE', 'ilvA', 'serA', 'serB', 'slyD']
    bitscore_per_pos = 0.7

    print(f"{'=' * 60}")
    print(f"Submitting MSA regeneration jobs")
    print(f"Bitscore-per-pos: {bitscore_per_pos}")
    print(f"Genes: {', '.join(genes)}")
    print(f"{'=' * 60}\n")

    # Submit all jobs
    job_map = {}  # {job_id: gene}
    for gene in genes:
        job_id = submit_msa_job(gene, bitscore_per_pos)
        if job_id:
            job_map[job_id] = gene

    if not job_map:
        print("\nERROR: No jobs were successfully submitted")
        sys.exit(1)

    print(f"\n{'=' * 60}")
    print(f"Submitted {len(job_map)} jobs")
    print(f"Job IDs: {', '.join(job_map.keys())}")
    print(f"{'=' * 60}\n")

    # Monitor jobs
    print("Monitoring job status (checking every 30 seconds)...")
    print("Press Ctrl+C to stop monitoring (jobs will continue running)\n")

    try:
        completed = set()
        while len(completed) < len(job_map):
            status_map = check_job_status(list(job_map.keys()))

            # Check which jobs are no longer in queue (completed or failed)
            for job_id, gene in job_map.items():
                if job_id not in status_map and job_id not in completed:
                    completed.add(job_id)
                    print(f"✓ Job {job_id} ({gene}) completed")

            # Show running/pending jobs
            running = [f"{job_id} ({job_map[job_id]}): {status}"
                      for job_id, status in status_map.items()]
            if running:
                print(f"[{time.strftime('%H:%M:%S')}] Active jobs: {', '.join(running)}")

            if len(completed) < len(job_map):
                time.sleep(30)

    except KeyboardInterrupt:
        print("\n\nMonitoring stopped. Jobs are still running in the background.")
        print("Check status with: squeue -u $USER")
        print("View logs in: logs/regen_msa_*.out")
        sys.exit(0)

    print(f"\n{'=' * 60}")
    print(f"All {len(job_map)} jobs completed!")
    print(f"{'=' * 60}")

    # Check which MSAs were successfully generated
    print("\nVerifying generated MSAs...")
    success_count = 0
    for gene in genes:
        aln_file = Path(f"data/msa/{gene}.aln")
        if aln_file.exists():
            with open(aln_file) as f:
                n_seqs = sum(1 for line in f if line.startswith('>'))
            print(f"✓ {gene}: {n_seqs} sequences")
            success_count += 1
        else:
            print(f"✗ {gene}: MSA file not found")

    print(f"\n{'=' * 60}")
    print(f"Summary: {success_count}/{len(genes)} MSAs generated successfully")
    print(f"{'=' * 60}")

    if success_count == len(genes):
        print("\nNext steps:")
        print("1. Transfer data/msa/*.aln to GPU cluster")
        print("2. Run ccmpred on each:")
        print("   for gene in metR cysB trpD trpE cysE ilvA serA serB slyD; do")
        print("     ccmpred -n 100 -d 0 data/msa/${gene}.aln data/raw/${gene}.raw")
        print("   done")
        print("3. Transfer .raw files back to data/raw/")
    else:
        print("\nCheck logs for failed jobs: logs/regen_msa_*.err")

    sys.exit(0 if success_count == len(genes) else 1)



if __name__ == "__main__":
    main()
