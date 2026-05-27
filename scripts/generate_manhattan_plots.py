"""
Generate Manhattan plots for the demo CSV outputs.

This script reads the sorted and thinned CSV outputs produced by the
Python and C++ pipelines and generates Manhattan plots for each dataset.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from thin_sorted_pvalues import logarithmic_thinning


matplotlib.use("Agg")


# Chromosome lengths taken from GRCh38; these are used to keep a consistent
# relative layout across raw and thinned plots.
CHROM_LENGTHS = {
    1: 248_956_422,
    2: 242_193_529,
    3: 198_295_559,
    4: 190_214_555,
    5: 181_538_259,
    6: 170_805_979,
    7: 159_345_973,
    8: 145_138_636,
    9: 138_394_717,
    10: 133_797_422,
    11: 135_086_622,
    12: 133_275_309,
    13: 114_364_328,
    14: 107_043_718,
    15: 101_991_189,
    16: 90_338_345,
    17: 83_257_441,
    18: 80_373_285,
    19: 58_617_616,
    20: 64_444_167,
    21: 46_709_983,
    22: 50_818_468,
}
CHROM_ORDER = tuple(sorted(CHROM_LENGTHS))


@dataclass(frozen=True)
class PlotSpec:
    input_csv: Path
    output_png: Path
    title: str
    requires_thinning: bool = False
    sorted_source: Path | None = None
    thinning_factor: float = 1.25
    annotate_count: bool = False


def canonical_sort(df: pd.DataFrame) -> pd.DataFrame:
    """Match the canonical ordering used throughout the regression tests."""
    df = df.copy()
    df["pp"] = df["pp"].round(6)
    df["original_index"] = df["original_index"].astype(int)
    df["chr"] = df["chr"].astype(int)
    df["pos"] = df["pos"].astype(int)
    return df.sort_values(
        by=["chr", "pos", "pp", "original_index"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)


def canonicalize_thinned(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure thinned outputs have consistent dtypes and order."""
    df = df.copy()
    df["pp"] = df["pp"].round(6)
    df["row_number"] = df["row_number"].astype(int)
    df["original_index"] = df["original_index"].astype(int)
    df["chr"] = df["chr"].astype(int)
    df["pos"] = df["pos"].astype(int)
    return df.sort_values("row_number").reset_index(drop=True)


def apply_thinning(sorted_df: pd.DataFrame, factor: float) -> pd.DataFrame:
    """Mirror the thinning logic from the Python/C++ pipeline."""
    working = sorted_df.copy()
    working["pp_key"] = (working["pp"] * 1_000_000).round().astype(np.int64)

    significance = working.sort_values(
        by=["pp_key", "chr", "pos", "original_index"],
        ascending=[False, True, True, True],
    ).reset_index()
    significance["rank"] = np.arange(1, len(significance) + 1, dtype=np.int64)

    rows_to_keep = logarithmic_thinning(len(significance), factor)
    rows_to_keep = [row for row in rows_to_keep if row <= len(significance)]
    top_n = min(len(significance), 2000)
    rows_to_keep = sorted(set(rows_to_keep).union(range(1, top_n + 1)))

    keep_mask = significance["rank"].isin(rows_to_keep)
    kept = significance.loc[keep_mask, ["index", "rank"]].sort_values("rank")

    thinned = sorted_df.iloc[kept["index"].to_numpy()].copy()
    thinned["row_number"] = kept["rank"].to_numpy()
    return canonicalize_thinned(thinned)


def _compute_positions(df: pd.DataFrame) -> Tuple[pd.Series, np.ndarray, Tuple[str, ...], float]:
    """
    Convert chromosome/position pairs into a consistent axis where each chromosome
    occupies a span proportional to its actual length in base pairs.
    """
    offsets = {}
    xticks = []
    tick_labels = []
    cumulative = 0
    for chrom in CHROM_ORDER:
        length = CHROM_LENGTHS[chrom]
        offsets[chrom] = cumulative
        xticks.append((cumulative + length / 2) / 1e6)
        tick_labels.append(str(chrom))
        cumulative += length

    if not set(df["chr"]).issubset(offsets):
        missing = set(df["chr"]) - set(offsets)
        raise ValueError(f"Encountered chromosomes without configured lengths: {sorted(missing)}")

    plot_positions = (df["pos"] + df["chr"].map(offsets)) / 1e6

    return plot_positions, np.array(xticks), tuple(tick_labels), cumulative / 1e6


def _load_dataframe(spec: PlotSpec) -> pd.DataFrame:
    if spec.requires_thinning:
        if spec.sorted_source is None:
            raise ValueError(f"{spec.title} requires a sorted source CSV to reapply thinning.")
        sorted_df = canonical_sort(pd.read_csv(spec.sorted_source))
        thinned_df = apply_thinning(sorted_df, spec.thinning_factor)
        # Persist regenerated thinning so downstream consumers stay in sync.
        spec.input_csv.parent.mkdir(parents=True, exist_ok=True)
        thinned_df.to_csv(spec.input_csv, index=False)
        return thinned_df

    return canonical_sort(pd.read_csv(spec.input_csv))


def _draw_manhattan(spec: PlotSpec) -> None:
    df = _load_dataframe(spec)
    if df.empty:
        raise ValueError(f"{spec.input_csv} is empty; cannot plot.")

    plot_positions, xticks, tick_labels, x_max = _compute_positions(df)
    colors = np.where(df["chr"] % 2 == 0, "#4C72B0", "#DD8452")

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(plot_positions, df["pp"], c=colors, s=6, linewidth=0, alpha=0.8)

    ax.set_title(spec.title)
    ax.set_xlabel("Genomic position (Mb)")
    ax.set_ylabel("log10 p-value")
    ax.set_xticks(xticks)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_xlim(0, x_max)
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    ax.margins(x=0.01)
    if spec.annotate_count:
        ax.text(
            0.99,
            0.96,
            f"n = {len(df):,}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"),
        )

    spec.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(spec.output_png, dpi=150)
    plt.close(fig)


def generate_plots(specs: Iterable[PlotSpec]) -> None:
    for spec in specs:
        _draw_manhattan(spec)


def _default_specs(repo_root: Path) -> Tuple[PlotSpec, ...]:
    demo_dir = repo_root / "demo_outputs"
    python_sorted = demo_dir / "python_final_sorted_data.csv"
    cpp_sorted = demo_dir / "cpp_final_sorted_data.csv"

    return (
        PlotSpec(
            input_csv=python_sorted,
            output_png=demo_dir / "python_raw_manhattan.png",
            title="Python Pipeline – Full Data",
            annotate_count=True,
        ),
        PlotSpec(
            input_csv=cpp_sorted,
            output_png=demo_dir / "cpp_raw_manhattan.png",
            title="C++ Pipeline – Full Data",
            annotate_count=True,
        ),
        PlotSpec(
            input_csv=demo_dir / "python_thinned.csv",
            output_png=demo_dir / "python_thinned_manhattan.png",
            title="Python Pipeline – Thinned",
            requires_thinning=False,
            sorted_source=None,
            annotate_count=True,
        ),
        PlotSpec(
            input_csv=demo_dir / "cpp_thinned.csv",
            output_png=demo_dir / "cpp_thinned_manhattan.png",
            title="C++ Pipeline – Thinned",
            requires_thinning=False,
            sorted_source=None,
            annotate_count=True,
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Manhattan plots for pipeline outputs."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Location of the repository root (defaults to scripts/../..).",
    )
    parser.add_argument(
        "--thinning-factor",
        type=float,
        default=1.0003,
        help="Logarithmic thinning factor to use when regenerating thinned datasets.",
    )
    args = parser.parse_args()

    specs = _default_specs(args.repo_root)
    # Allow overriding the thinning factor without redefining specs manually.
    specs = tuple(
        spec if not spec.requires_thinning else PlotSpec(
            input_csv=spec.input_csv,
            output_png=spec.output_png,
            title=spec.title,
            requires_thinning=True,
            sorted_source=spec.sorted_source,
            thinning_factor=args.thinning_factor,
            annotate_count=spec.annotate_count,
        )
        for spec in specs
    )
    generate_plots(specs)


if __name__ == "__main__":
    main()
