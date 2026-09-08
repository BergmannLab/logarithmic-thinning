#!/usr/bin/env python3
"""
Variant of thin_sorted_pvalues that also writes the retained indices to a CSV.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Iterable

from thin_sorted_pvalues import (
    DEFAULT_METHOD,
    DEFAULT_THINNING_FACTOR,
    IO_BUFFER_SIZE,
    UNTHINNED_COLUMN,
    filter_specific_rows_awk,
    filter_specific_rows_python,
    logarithmic_thinning,
    time_execution,
)

ROWS_CSV_LABEL = "Row_Index"


def save_rows_to_csv(rows_to_keep: Iterable[int], output_csv: str) -> None:
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([ROWS_CSV_LABEL])
        for row in rows_to_keep:
            writer.writerow([row])
    print(f"Rows to keep saved to {output_csv}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Logarithmically thin a sorted file and emit retained indices."
    )
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
    parser.add_argument(
        "--rows-output",
        help="Path to write the retained row indices CSV. "
        "Defaults to rows_to_keep_{input}.csv in the same directory.",
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

    base_dir = os.path.dirname(input_file)

    if args.output:
        output_file = args.output
    else:
        suffix = method.lower()
        output_file = os.path.join(
            base_dir, f"thinned_{suffix}_{os.path.basename(input_file)}"
        )

    if args.rows_output:
        rows_csv_file = args.rows_output
    else:
        rows_csv_file = os.path.join(
            base_dir, f"rows_to_keep_{os.path.basename(input_file)}.csv"
        )

    with open(input_file, "r", encoding="utf-8") as handle:
        total_rows = sum(1 for _ in handle)
    if total_rows <= 1:
        raise ValueError(f"Input file {input_file} does not contain data rows.")
    data_rows = total_rows - 1
    print(f"Detected {data_rows} data rows (plus header) in {input_file}")

    rows_to_keep = logarithmic_thinning(data_rows, thinning_factor)
    save_rows_to_csv(rows_to_keep, rows_csv_file)

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
