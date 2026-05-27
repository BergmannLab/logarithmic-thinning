#!/dcsrsoft/spack/hetre/v1.2/spack/opt/spack/linux-centos7-sandybridge/gcc-9.3.0/python-3.8.8-7p6ikdiracpfzkqxrpo6i274b7dxe2ga/bin/python
#SBATCH --account=<account>
#SBATCH --job-name=iterative_merge_sort
#SBATCH --error=slurm-%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=2:00:00
#SBATCH --partition=urblauna

import os
import numpy as np
import pandas as pd
from heapq import merge
from multiprocessing import Pool
import time

def merge_pair(args):
    """Merge a pair of sorted files using heapq.merge."""
    file1, file2, output_dir, idx = args

    print(f"Merging {file1} and {file2}")

    # Load the data from both files
    data1 = np.load(file1)
    data2 = np.load(file2)

    pp1, indices1, chr1, pos1 = data1['pp'], data1['indices'], data1['chr'], data1['pos']
    pp2, indices2, chr2, pos2 = data2['pp'], data2['indices'], data2['chr'], data2['pos']

    # Merge the data using heapq.merge (already sorted input)
    merged = list(merge(
        zip(-pp1, indices1, chr1, pos1),
        zip(-pp2, indices2, chr2, pos2)
    ))

    merged_pp, merged_indices, merged_chr, merged_pos = zip(*[(-pp, idx, ch, pos) for pp, idx, ch, pos in merged])

    # Convert to numpy arrays
    merged_pp = np.array(merged_pp, dtype=np.float32)
    merged_indices = np.array(merged_indices, dtype=np.int32)
    merged_chr = np.array(merged_chr, dtype=np.int8)
    merged_pos = np.array(merged_pos, dtype=np.int32)

    # Save the merged chunk
    merged_file = os.path.join(output_dir, f"merged_{time.time_ns()}.npz")
    np.savez_compressed(merged_file, pp=merged_pp, indices=merged_indices, chr=merged_chr, pos=merged_pos)

    # Clean up
    os.remove(file1)
    os.remove(file2)

    return merged_file

def merge_chunks(input_files, output_dir):
    """
    Iteratively merge sorted chunks by pairing the largest with the smallest file using heapq.merge
    in parallel until a single sorted list remains. Save the final list as both a binary numpy file and a CSV.
    """
    while len(input_files) > 1:
        new_files = []

        # Pair files: largest with smallest
        input_files = sorted(input_files, key=lambda f: os.path.getsize(f))
        pairs = [(input_files[i], input_files[-(i + 1)], output_dir, i) for i in range(len(input_files) // 2)]

        # Use multiprocessing to merge pairs in parallel
        with Pool() as pool:
            new_files.extend(pool.map(merge_pair, pairs))

        # If there's an odd number of files, carry the middle one forward
        if len(input_files) % 2 == 1:
            new_files.append(input_files[len(input_files) // 2])

        input_files = new_files

    # Final file
    final_file = input_files[0]
    print(f"Final sorted file: {final_file}")

    # Load final file and save in desired formats
    final_data = np.load(final_file)
    final_pp = final_data['pp']
    final_indices = final_data['indices']
    final_chr = final_data['chr']
    final_pos = final_data['pos']

    # Save as binary numpy format
    binary_file = os.path.join(output_dir, "final_sorted_data.npy")
    np.save(binary_file, {
        'pp': final_pp,
        'indices': final_indices,
        'chr': final_chr,
        'pos': final_pos
    })

    # Save as CSV
    csv_file = os.path.join(output_dir, "final_sorted_data.csv")
    df = pd.DataFrame({
        'chr': final_chr,
        'pos': final_pos,
        'pp': final_pp,
        'original_index': final_indices
    })
    df.to_csv(csv_file, index=False)

    print(f"Saved final data to {binary_file} and {csv_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Iteratively merge sorted chunks.")
    parser.add_argument("--input_dir", required=True, help="Directory containing sorted chunk files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the final sorted output.")

    args = parser.parse_args()

    # Get sorted chunk files
    input_files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if f.endswith(".npz")]

    # Perform iterative merge
    merge_chunks(input_files, args.output_dir)

