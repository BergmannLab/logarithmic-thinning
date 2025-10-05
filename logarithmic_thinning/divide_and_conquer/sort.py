#!/usr/bin/env python3
import os
import sys
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
import time

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process chunks and output sorted data.')
    parser.add_argument('input_file', help='Path to the input file.')
    parser.add_argument('--chunksize', type=int, default=10000000000, help='Number of rows per chunk.')
    parser.add_argument('--shared_dir', required=True, help='Path to the shared directory for output chunks.')
    parser.add_argument('--pp_threshold', required=True, help='Smallest -log10 p-value to keep.')
    return parser.parse_args()

def get_log10p_columns(input_file):
    with open(input_file, 'r') as f:
        headers = f.readline().strip().split()
    log10p_columns = [col for col in headers if col.endswith('true-log10p')]
    print(f"Detected log10p columns: {log10p_columns}", flush=True)
    return log10p_columns

def process_chunk(args):
    chunk_data, start_row, shared_dir, log10p_columns, pp_threshold = args
    n_rows, n_cols = chunk_data[log10p_columns].shape
    total_elements = n_rows * n_cols  # Store the total size of the data

    print(f"Processing chunk: {n_rows} rows, {n_cols} cols, total elements: {total_elements}", flush=True)

    # Flatten `pp` values and generate associated data
    combined_pp = chunk_data[log10p_columns].values.flatten(order='F').astype(np.float32)
    row_indices = np.tile(np.arange(start_row, start_row + n_rows, dtype=np.int32), n_cols)
    chr_values = np.repeat(chunk_data['chr'].values.astype(np.int8), n_cols)
    pos_values = np.repeat(chunk_data['pos'].values.astype(np.int32), n_cols)
    
    # Delate chunk data, making enabling memory release
    del chunk_data

    # Replace NaNs
    combined_pp = np.nan_to_num(combined_pp, nan=-np.inf)

    # Apply threshold filter
    mask = combined_pp > pp_threshold
    filtered_pp = combined_pp[mask]
    filtered_row_indices = row_indices[mask]
    filtered_chr = chr_values[mask]
    filtered_pos = pos_values[mask]

    # Sort filtered data
    sorted_indices = np.argsort(-filtered_pp)
    sorted_pp = filtered_pp[sorted_indices]
    sorted_row_indices = filtered_row_indices[sorted_indices].astype(np.int32)
    sorted_chr = filtered_chr[sorted_indices]
    sorted_pos = filtered_pos[sorted_indices]

    print(len(mask), mask)

    # Verify the size of the sorted data
    assert len(sorted_pp) == len(sorted_row_indices) == len(sorted_chr) == len(sorted_pos), "Mismatch in data size after filtering and sorting"

    # Save sorted data to a unique file in the shared directory
    unique_id = f"{time.time_ns()}_{os.getpid()}"
    chunk_filename = os.path.join(shared_dir, f"chunk_{unique_id}.npz")
    np.savez_compressed(chunk_filename, pp=sorted_pp, indices=sorted_row_indices, chr=sorted_chr, pos=sorted_pos)
    print(f"Saved chunk to {chunk_filename}", flush=True)

    return chunk_filename

def main():
    args = parse_arguments()
    input_file = args.input_file
    chunksize = args.chunksize
    shared_dir = args.shared_dir
    pp_threshold = args.pp_threshold
    print(pp_threshold, flush=True)
    pp_threshold = int(pp_threshold)
    print(pp_threshold, flush=True)
    if not os.path.exists(shared_dir):
        os.makedirs(shared_dir)

    num_cpus = mp.cpu_count()
    print(f"Using {num_cpus} CPUs.", flush=True)

    print(f"Input file: {input_file}.", flush=True)

    log10p_columns = get_log10p_columns(input_file) 
    if not log10p_columns:
        print('No columns ending with log10p found.', flush=True)
        sys.exit(1)

    all_columns = ['chr', 'pos'] + log10p_columns
    dtype_dict = {col: np.float32 for col in log10p_columns}
    dtype_dict.update({'chr': np.int8, 'pos': np.int32})

    reader = pd.read_csv(
        input_file,
        sep=' ',
        usecols=all_columns,
        chunksize=chunksize,
        dtype=dtype_dict,
        na_values=['NA', 'NaN', '']
    )

    pool = mp.Pool(processes=num_cpus)
    start_row = 0
    tasks = []

    for chunk in reader:
        print(f"Processing chunk starting at row {start_row}", flush=True)
        task_args = (chunk, start_row, shared_dir, log10p_columns, pp_threshold)
        task = pool.apply_async(process_chunk, args=(task_args,))
        tasks.append(task)
        start_row += len(chunk)

    pool.close()
    pool.join()

    print("All chunks have been processed and saved.", flush=True)

if __name__ == '__main__':
    main()

