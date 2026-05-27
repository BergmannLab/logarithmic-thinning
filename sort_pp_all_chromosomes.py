#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import time
import heapq

def parse_arguments():
    parser = argparse.ArgumentParser(description='Wrapper script to run maxisort.py on all chromosomes and merge results.')
    parser.add_argument('input_folder', help='Path to the input folder containing chromosome files.')
    parser.add_argument('output_prefix', help='Prefix for the final merged output files.')
    parser.add_argument('--max_cpus', type=int, default=47, help='Maximum number of CPUs to allocate for each job.')
    parser.add_argument('--chunksize', type=int, default=100000, help='Chunksize for maxisort.py.')
    parser.add_argument('--wait_time', type=int, default=30, help='Time to wait between job status checks (in seconds).')
    return parser.parse_args()

def submit_maxisort_job(chromosome, input_file, output_prefix, max_cpus, chunksize):
    job_name = f"maxisort_chr{chromosome}"
    slurm_script = f"""#!/bin/bash
#SBATCH --account=<account>
#SBATCH --job-name={job_name}
#SBATCH --output=slurm-{job_name}-%j.out
#SBATCH --error=slurm-{job_name}-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={max_cpus}
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --partition=urblauna

python -u maxisort.py {input_file} /{output_prefix}_chr{chromosome} --chunksize {chunksize}
"""
    slurm_filename = f"slurm_maxisort_chr{chromosome}.sh"
    with open(slurm_filename, 'w') as f:
        f.write(slurm_script)
    # Submit the job
    result = subprocess.run(['sbatch', slurm_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted job for chromosome {chromosome}: Job ID {job_id}", flush=True)
        return job_id
    else:
        print(f"Failed to submit job for chromosome {chromosome}: {result.stderr}", flush=True)
        sys.exit(1)

def check_job_status(job_ids):
    completed_jobs = set()
    # Use sacct to check job status
    result = subprocess.run(['sacct', '-n', '-X', '--format=JobID,State', '--jobs', ','.join(job_ids)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        lines = result.stdout.strip().split('\n')
        for line in lines:
            job_id, state = line.strip().split()
            if state in ['COMPLETED']:
                completed_jobs.add(job_id.split('.')[0])
            elif state in ['FAILED', 'CANCELLED', 'TIMEOUT', 'NODE_FAIL', 'PREEMPTED']:
                print(f"Job {job_id} failed with state {state}", flush=True)
                sys.exit(1)
    else:
        print(f"Failed to check job status: {result.stderr}", flush=True)
        sys.exit(1)
    return completed_jobs

def merge_files(file1_prefix, file2_prefix, output_prefix):
    # Merging two sorted pp and index files using heapq
    pp_file1 = f"{file1_prefix}_pp.npy"
    idx_file1 = f"{file1_prefix}_original_python_index.npy"
    pp_file2 = f"{file2_prefix}_pp.npy"
    idx_file2 = f"{file2_prefix}_original_python_index.npy"

    output_pp_file = f"{output_prefix}_pp.npy"
    output_idx_file = f"{output_prefix}_original_python_index.npy"

    # Load numpy arrays
    import numpy as np
    pp1 = np.load(pp_file1, mmap_mode='r')
    idx1 = np.load(idx_file1, mmap_mode='r')
    pp2 = np.load(pp_file2, mmap_mode='r')
    idx2 = np.load(idx_file2, mmap_mode='r')

    # Initialize heap
    heap = []
    i, j = 0, 0
    if len(pp1) > 0:
        heapq.heappush(heap, (-pp1[0], 0, 'pp1'))
    if len(pp2) > 0:
        heapq.heappush(heap, (-pp2[0], 0, 'pp2'))

    merged_pp = []
    merged_idx = []

    while heap:
        neg_pp_val, idx, source = heapq.heappop(heap)
        pp_val = -neg_pp_val
        if source == 'pp1':
            merged_pp.append(pp_val)
            merged_idx.append(idx1[idx])
            i = idx + 1
            if i < len(pp1):
                heapq.heappush(heap, (-pp1[i], i, 'pp1'))
        else:
            merged_pp.append(pp_val)
            merged_idx.append(idx2[idx])
            j = idx + 1
            if j < len(pp2):
                heapq.heappush(heap, (-pp2[j], j, 'pp2'))

    # Save merged arrays
    np.save(output_pp_file, np.array(merged_pp, dtype=np.float64))
    np.save(output_idx_file, np.array(merged_idx, dtype=np.int64))

    # Clean up individual files
    os.remove(pp_file1)
    os.remove(idx_file1)
    os.remove(pp_file2)
    os.remove(idx_file2)

def main():
    args = parse_arguments()
    input_folder = args.input_folder
    output_prefix = args.output_prefix
    max_cpus = args.max_cpus
    chunksize = args.chunksize
    wait_time = args.wait_time

    output_prefix = input_folder+"/"+output_prefix

    # Submit jobs for all chromosomes
    job_ids = {}
    for chrom in range(1, 23):
        input_file = os.path.join(input_folder, f"output_ukb_imp_chr{chrom}_v3.txt")
        if not os.path.exists(input_file):
            print(f"Input file for chromosome {chrom} does not exist: {input_file}", flush=True)
            sys.exit(1)
        job_id = submit_maxisort_job(chrom, input_file, output_prefix, max_cpus, chunksize)
        job_ids[str(job_id)] = chrom

    # Monitor job completion
    completed_chromosomes = []
    while len(completed_chromosomes) < 22:
        time.sleep(wait_time)
        completed_jobs = check_job_status(list(job_ids.keys()))
        for job_id in completed_jobs:
            chrom = job_ids[job_id]
            if chrom not in completed_chromosomes:
                print(f"Chromosome {chrom} processing completed.", flush=True)
                completed_chromosomes.append(chrom)
                # As soon as at least two chromosomes are processed, start merging
                if len(completed_chromosomes) >= 2:
                    # Merge available chromosomes
                    available_chroms = sorted(completed_chromosomes)
                    while len(available_chroms) >= 2:
                        chrom1 = available_chroms.pop(0)
                        chrom2 = available_chroms.pop(0)
                        print(f"Merging chromosomes {chrom1} and {chrom2}", flush=True)
                        prefix1 = f"{output_prefix}_chr{chrom1}"
                        prefix2 = f"{output_prefix}_chr{chrom2}"
                        merged_prefix = f"{output_prefix}_merged_chr{chrom1}_{chrom2}"
                        merge_files(prefix1, prefix2, merged_prefix)
                        # Add the merged file back to the list for further merging
                        merged_chrom = f"merged_chr{chrom1}_{chrom2}"
                        completed_chromosomes.append(merged_chrom)
                        available_chroms.append(merged_chrom)
        # Remove processed jobs from job_ids
        for job_id in completed_jobs:
            if job_id in job_ids:
                del job_ids[job_id]

    # Final merging if necessary
    if len(completed_chromosomes) == 1:
        final_prefix = f"{output_prefix}_pp.npy"
        os.rename(f"{output_prefix}_{completed_chromosomes[0]}_pp.npy", final_prefix)
        idx_prefix = f"{output_prefix}_original_python_index.npy"
        os.rename(f"{output_prefix}_{completed_chromosomes[0]}_original_python_index.npy", idx_prefix)
        print(f"Final sorted pp and index files are saved with prefix {output_prefix}", flush=True)
    else:
        # Merge remaining files
        while len(completed_chromosomes) > 1:
            chrom1 = completed_chromosomes.pop(0)
            chrom2 = completed_chromosomes.pop(0)
            print(f"Merging chromosomes {chrom1} and {chrom2}", flush=True)
            prefix1 = f"{output_prefix}_{chrom1}"
            prefix2 = f"{output_prefix}_{chrom2}"
            merged_prefix = f"{output_prefix}_merged_{chrom1}_{chrom2}"
            merge_files(prefix1, prefix2, merged_prefix)
            merged_chrom = f"merged_{chrom1}_{chrom2}"
            completed_chromosomes.append(merged_chrom)
        # Rename final merged files
        final_prefix = f"{output_prefix}_pp.npy"
        os.rename(f"{output_prefix}_{completed_chromosomes[0]}_pp.npy", final_prefix)
        idx_prefix = f"{output_prefix}_original_python_index.npy"
        os.rename(f"{output_prefix}_{completed_chromosomes[0]}_original_python_index.npy", idx_prefix)
        print(f"Final sorted pp and index files are saved with prefix {output_prefix}", flush=True)

if __name__ == "__main__":
    main()

