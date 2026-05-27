#!/bin/python3
#SBATCH --account=<account>
#SBATCH --job-name=thin
#SBATCH --error=slurm-%j.err
#SBATCH --output=slurm-%j.out
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 100G
#SBATCH --time 1:00:00
#SBATCH --partition urblauna

"""
Utility to apply logarithmic thinning to a pre-sorted table.
Configuration defaults live at the top of the file for quick inspection.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import tempfile
import time
from typing import Iterable

# ---------------------------------------------------------------------------
# User-facing defaults
# ---------------------------------------------------------------------------

DEFAULT_THINNING_FACTOR = 1.0003
DEFAULT_METHOD = "python"
UNTHINNED_COLUMN = "unthinned_rank"
IO_BUFFER_SIZE = 16 * 1024 * 1024


# ---------------------------------------------------------------------------
# Core thinning helpers
# ---------------------------------------------------------------------------

def determine_threshold(factor: float) -> int:
    return round(factor / (factor - 1.0))


def logarithmic_thinning(total_rows: int, factor: float) -> list[int]:
    print("Number of p-values before thinning:", total_rows)
    indices: list[int] = []
    index = total_rows
    indices.append(index)

    threshold = determine_threshold(factor)
    print("Number of unthinned top p-values:", threshold)

    while index > threshold:
        index = math.floor(index / factor)
        indices.append(index)

    indices.extend(range(max(threshold - 1, 1), 0, -1))
    indices = sorted(set(indices))
    print("Total number of p-values after thinning:", len(indices))
    return indices


def filter_specific_rows_python(
    input_file: str, output_file: str, rows_to_keep: Iterable[int]
) -> None:
    rows_set = set(rows_to_keep)
    with open(input_file, "r", buffering=IO_BUFFER_SIZE) as infile, open(
        output_file, "w", buffering=IO_BUFFER_SIZE
    ) as outfile:
        header = infile.readline().rstrip("\n")
        outfile.write(f"{header},{UNTHINNED_COLUMN}\n")

        for data_idx, line in enumerate(infile, start=1):
            if data_idx in rows_set:
                outfile.write(line.rstrip("\n") + f",{data_idx}\n")


def filter_specific_rows_awk(
    input_file: str, output_file: str, rows_to_keep: Iterable[int]
) -> None:
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
        for row in rows_to_keep:
            temp_file.write(f"{row}\n")
        temp_filename = temp_file.name

    awk_command = (
        "awk 'NR==FNR { lines[$1]; next } "
        f"FNR==1 {{ print $0\",{UNTHINNED_COLUMN}\"; next }} "
        "(FNR-1) in lines { print $0\",\"(FNR-1) }' "
        f"{temp_filename} {input_file} > {output_file}"
    )
    subprocess.run(awk_command, shell=True, check=True)
    os.remove(temp_filename)


def time_execution(func, *args) -> float:
    start = time.time()
    func(*args)
    end = time.time()
    return end - start


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Apply logarithmic thinning to a file.")
    parser.add_argument("input_file", help="Path to the input (sorted) file.")
    parser.add_argument(
        "pos_thinning_factor",
        nargs="?",
        type=float,
        help="(Deprecated positional) thinning factor.",
    )
    parser.add_argument(
        "pos_method",
        nargs="?",
        choices=["python", "awk"],
        help="(Deprecated positional) method.",
    )
    parser.add_argument(
        "--thinning-factor",
        dest="opt_thinning_factor",
        type=float,
        default=DEFAULT_THINNING_FACTOR,
        help=f"Logarithmic thinning factor (default: {DEFAULT_THINNING_FACTOR}).",
    )
    parser.add_argument(
        "--method",
        dest="opt_method",
        choices=["python", "awk"],
        default=DEFAULT_METHOD,
        help=f"Row selection method (default: {DEFAULT_METHOD}).",
    )
    parser.add_argument(
        "--output",
        help="Explicit output path. Defaults to thinned_{method}_{input}.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    input_file = args.input_file
    thinning_factor = (
        args.pos_thinning_factor
        if args.pos_thinning_factor is not None
        else args.opt_thinning_factor
    )
    method = (
        args.pos_method if args.pos_method is not None else args.opt_method
    )

    if thinning_factor <= 1.0:
        raise ValueError("Thinning factor must be greater than 1.0")

    if args.output:
        output_file = args.output
    else:
        suffix = method.lower()
        output_dir = os.path.dirname(input_file)
        output_name = f"thinned_{suffix}_{os.path.basename(input_file)}"
        output_file = os.path.join(output_dir, output_name)

    with open(input_file, "r", encoding="utf-8") as handle:
        total_rows = sum(1 for _ in handle)
    if total_rows <= 1:
        raise ValueError(f"Input file {input_file} does not contain data rows.")
    data_rows = total_rows - 1
    print(f"Detected {data_rows} data rows (plus header) in {input_file}")

    rows_to_keep = logarithmic_thinning(data_rows, thinning_factor)

    if method == "python":
        print("Running Python version...")
        elapsed = time_execution(
            filter_specific_rows_python, input_file, output_file, rows_to_keep
        )
        print(f"Python version took {elapsed:.2f} seconds.")
    else:
        print("Running AWK version...")
        elapsed = time_execution(
            filter_specific_rows_awk, input_file, output_file, rows_to_keep
        )
        print(f"AWK version took {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
