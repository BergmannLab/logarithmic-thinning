#!/usr/bin/env python3
"""
Chunked sorter for GWAS summary statistics (-log10 p-values).

Key configuration knobs operate at the top of this file so they are easy to audit.
All options remain overridable from the command line.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from typing import Iterable, List, Sequence

import multiprocessing as mp
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# User-facing defaults (can be overridden via CLI)
# ---------------------------------------------------------------------------

DEFAULT_CHUNK_SIZE = 5_000_000
DEFAULT_PP_THRESHOLD = 0.0  # keep everything unless specified otherwise
DEFAULT_DELIMITER = " "
DEFAULT_CHR_COLUMN = "chr"
DEFAULT_POS_COLUMN = "pos"
DEFAULT_PP_COLUMNS: Sequence[str] = ()
DEFAULT_P_COLUMNS: Sequence[str] = ()


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def _split_columns(value: str | None) -> List[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Process chunks and output sorted data.")
    parser.add_argument("input_file", help="Path to the input summary statistics file.")
    parser.add_argument("--shared_dir", required=True, help="Directory for chunk outputs.")
    parser.add_argument(
        "--chunksize",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help=f"Number of rows per chunk (default: {DEFAULT_CHUNK_SIZE:,}).",
    )
    parser.add_argument(
        "--delimiter",
        default=DEFAULT_DELIMITER,
        help=f"Field delimiter in the input file (default: '{DEFAULT_DELIMITER}').",
    )
    parser.add_argument(
        "--chr-column",
        default=DEFAULT_CHR_COLUMN,
        help=f"Chromosome column name (default: {DEFAULT_CHR_COLUMN!r}).",
    )
    parser.add_argument(
        "--pos-column",
        default=DEFAULT_POS_COLUMN,
        help=f"Position column name (default: {DEFAULT_POS_COLUMN!r}).",
    )
    parser.add_argument(
        "--pp-columns",
        default=",".join(DEFAULT_PP_COLUMNS),
        help="Comma-separated list of columns already on the -log10(p) scale.",
    )
    parser.add_argument(
        "--p-columns",
        default=",".join(DEFAULT_P_COLUMNS),
        help="Comma-separated list of raw p-value columns (converted to -log10(p)).",
    )
    parser.add_argument(
        "--pp-threshold",
        type=float,
        default=DEFAULT_PP_THRESHOLD,
        help="Smallest -log10(p) value to retain (default: 0).",
    )
    parser.add_argument(
        "--auto-detect",
        action="store_true",
        help="Auto-detect columns ending with 'true-log10p' when --pp-columns is omitted.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SortConfig:
    shared_dir: str
    chr_column: str
    pos_column: str
    pp_columns: Sequence[str]
    p_columns: Sequence[str]
    threshold: float


def auto_detect_pp_columns(path: str, delimiter: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as handle:
        headers = handle.readline().strip().split(delimiter)
    detected = [col for col in headers if col.endswith("true-log10p")]
    print(f"Detected log10p columns: {detected}", flush=True)
    return detected


def process_chunk(args: tuple[pd.DataFrame, int, SortConfig]) -> str:
    chunk_data, start_row, config = args

    pp_mats: List[np.ndarray] = []

    if config.pp_columns:
        pp_mats.append(
            chunk_data[list(config.pp_columns)].to_numpy(dtype=np.float32, copy=True)
        )

    if config.p_columns:
        raw = chunk_data[list(config.p_columns)].to_numpy(dtype=np.float64, copy=True)
        raw = np.clip(raw, np.finfo(np.float64).tiny, 1.0)
        pp_mats.append((-np.log10(raw)).astype(np.float32))

    if not pp_mats:
        raise ValueError("No p-value columns supplied; nothing to process.")

    pp_matrix = np.concatenate(pp_mats, axis=1)
    n_rows, n_cols = pp_matrix.shape

    total_elements = n_rows * n_cols
    print(
        f"Processing chunk: {n_rows} rows × {n_cols} columns "
        f"(total {total_elements:,} elements)",
        flush=True,
    )

    combined_pp = pp_matrix.flatten(order="F")

    row_indices = np.tile(
        np.arange(start_row, start_row + n_rows, dtype=np.int32), n_cols
    )
    chr_values = np.tile(
        chunk_data[config.chr_column].to_numpy(dtype=np.int16, copy=False), n_cols
    )
    pos_values = np.tile(
        chunk_data[config.pos_column].to_numpy(dtype=np.int32, copy=False), n_cols
    )

    del chunk_data

    combined_pp = np.nan_to_num(combined_pp, nan=-np.inf)

    mask = combined_pp > config.threshold
    filtered_pp = combined_pp[mask]
    filtered_row_indices = row_indices[mask]
    filtered_chr = chr_values[mask]
    filtered_pos = pos_values[mask]

    if filtered_pp.size == 0:
        print("Chunk produced no values above threshold; skipping write.", flush=True)
        return ""

    sorted_indices = np.argsort(-filtered_pp)
    sorted_pp = filtered_pp[sorted_indices]
    sorted_row_indices = filtered_row_indices[sorted_indices].astype(np.int32)
    sorted_chr = filtered_chr[sorted_indices]
    sorted_pos = filtered_pos[sorted_indices]

    assert (
        len(sorted_pp)
        == len(sorted_row_indices)
        == len(sorted_chr)
        == len(sorted_pos)
    ), "Mismatch in data size after filtering and sorting"

    unique_id = f"{time.time_ns()}_{os.getpid()}"
    chunk_filename = os.path.join(config.shared_dir, f"chunk_{unique_id}.npz")
    np.savez_compressed(
        chunk_filename,
        pp=sorted_pp.astype(np.float32),
        indices=sorted_row_indices,
        chr=sorted_chr.astype(np.int16),
        pos=sorted_pos.astype(np.int32),
    )
    print(f"Saved chunk to {chunk_filename}", flush=True)
    return chunk_filename


def main() -> None:
    args = parse_arguments()

    pp_columns = _split_columns(args.pp_columns)
    p_columns = _split_columns(args.p_columns)

    if not pp_columns and args.auto_detect:
        pp_columns = auto_detect_pp_columns(args.input_file, args.delimiter)
    if not pp_columns and not p_columns:
        print(
            "ERROR: No p-value columns supplied. Use --pp-columns/--p-columns "
            "or enable --auto-detect.",
            flush=True,
        )
        sys.exit(1)

    os.makedirs(args.shared_dir, exist_ok=True)
    num_cpus = mp.cpu_count()
    print(f"Using {num_cpus} CPUs.", flush=True)
    print(f"Input file: {args.input_file}", flush=True)

    requested_columns: List[str] = [
        args.chr_column,
        args.pos_column,
        *pp_columns,
        *p_columns,
    ]

    dtype_dict = {
        **{col: np.float32 for col in pp_columns},
        **{col: np.float64 for col in p_columns},
        args.chr_column: np.int32,
        args.pos_column: np.int32,
    }

    reader = pd.read_csv(
        args.input_file,
        sep=args.delimiter,
        usecols=requested_columns,
        chunksize=args.chunksize,
        dtype=dtype_dict,
        na_values=["NA", "NaN", ""],
    )

    config = SortConfig(
        shared_dir=args.shared_dir,
        chr_column=args.chr_column,
        pos_column=args.pos_column,
        pp_columns=pp_columns,
        p_columns=p_columns,
        threshold=float(args.pp_threshold),
    )

    pool = mp.Pool(processes=num_cpus)
    start_row = 0

    for chunk in reader:
        print(f"Processing chunk starting at row {start_row}", flush=True)
        pool.apply_async(process_chunk, args=((chunk, start_row, config),))
        start_row += len(chunk)

    pool.close()
    pool.join()

    print("All chunks have been processed and saved.", flush=True)


if __name__ == "__main__":
    main()
