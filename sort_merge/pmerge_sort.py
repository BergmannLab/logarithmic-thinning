#!/usr/bin/env python3
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
import sys
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed, BrokenExecutor


def merge_pair(args):
    """Merge a pair of sorted .npz files using numpy concat + argsort."""
    file1, file2, output_dir, _ = args

    print(f"Merging {os.path.basename(file1)} + {os.path.basename(file2)}", flush=True)

    data1 = np.load(file1)
    data2 = np.load(file2)

    pp      = np.concatenate([data1["pp"],      data2["pp"]])
    indices = np.concatenate([data1["indices"],  data2["indices"]])
    chr_    = np.concatenate([data1["chr"],      data2["chr"]])
    pos     = np.concatenate([data1["pos"],      data2["pos"]])
    del data1, data2

    order = np.argsort(-pp, kind="stable")

    merged_file = os.path.join(output_dir, f"merged_{time.time_ns()}.npz")
    np.savez_compressed(
        merged_file,
        pp      = pp[order].astype(np.float32),
        indices = indices[order].astype(np.int32),
        chr     = chr_[order].astype(np.int16),
        pos     = pos[order].astype(np.int32),
    )
    del pp, indices, chr_, pos, order

    os.remove(file1)
    os.remove(file2)

    return merged_file


def merge_chunks(input_files, output_dir, n_workers):
    """
    Iteratively merge sorted chunks by pairing the largest with the smallest
    until a single sorted file remains.
    """
    round_n = 0
    while len(input_files) > 1:
        round_n += 1
        input_files = sorted(input_files, key=lambda f: os.path.getsize(f))
        n_pairs = len(input_files) // 2
        pairs = [
            (input_files[i], input_files[-(i + 1)], output_dir, i)
            for i in range(n_pairs)
        ]
        print(
            f"[round {round_n}] {len(input_files)} files -> {n_pairs} merges "
            f"(+{len(input_files) % 2} passthrough), workers={n_workers}",
            flush=True,
        )

        new_files = []
        try:
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                futs = {ex.submit(merge_pair, p): p for p in pairs}
                for fut in as_completed(futs):
                    new_files.append(fut.result())
        except BrokenExecutor as exc:
            print(
                f"ERROR: a merge worker died unexpectedly (likely OOM-killed). "
                f"({type(exc).__name__}: {exc})",
                file=sys.stderr, flush=True,
            )
            raise

        if len(input_files) % 2 == 1:
            new_files.append(input_files[n_pairs])

        input_files = new_files

    final_file = input_files[0]
    print(f"Final sorted file: {final_file}", flush=True)

    final_data = np.load(final_file)
    final_pp      = final_data["pp"]
    final_indices = final_data["indices"]
    final_chr     = final_data["chr"]
    final_pos     = final_data["pos"]

    binary_file = os.path.join(output_dir, "final_sorted_data.npy")
    np.save(binary_file, {
        "pp": final_pp, "indices": final_indices,
        "chr": final_chr, "pos": final_pos,
    })

    csv_file = os.path.join(output_dir, "final_sorted_data.csv")
    df = pd.DataFrame({
        "chr": final_chr, "pos": final_pos,
        "pp": final_pp, "original_index": final_indices,
    })
    df.to_csv(csv_file, index=False)
    print(f"Saved final data to {binary_file} and {csv_file}", flush=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Iteratively merge sorted chunks.")
    parser.add_argument("--input_dir",  required=True, help="Directory containing sorted chunk .npz files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the final sorted output.")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel merge workers. Defaults to $SLURM_CPUS_PER_TASK or 4.",
    )
    args = parser.parse_args()

    if args.workers is not None:
        n_workers = max(1, args.workers)
    else:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        n_workers = int(slurm_cpus) if slurm_cpus else 4

    input_files = [
        os.path.join(args.input_dir, f)
        for f in os.listdir(args.input_dir)
        if f.endswith(".npz")
    ]
    print(f"Found {len(input_files)} chunk files to merge, using {n_workers} workers.", flush=True)

    merge_chunks(input_files, args.output_dir, n_workers)
