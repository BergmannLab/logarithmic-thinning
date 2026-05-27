import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from tests.cpp_utils import REPO_ROOT, ensure_cpp_binaries

from thin_sorted_pvalues import logarithmic_thinning  # noqa: E402

BGENIE_HEADER = "chr pos traitA_true-log10p traitB_true-log10p\n"

# Chromosome lengths in base pairs (GRCh38 primary assembly)
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

CHROMOSOMES = tuple(CHROM_LENGTHS.keys())

BACKGROUND_SNPS = 100_000
THINNING_FACTOR = 1.002 # DEFAULT is 1.0003

PEAKS = [
    {"chrom": 3, "center_frac": 0.35, "height": 60.0, "decay": 1.2e6, "count": 140},
    {"chrom": 7, "center_frac": 0.55, "height": 45.0, "decay": 1.0e6, "count": 120},
    {"chrom": 19, "center_frac": 0.42, "height": 35.0, "decay": 0.7e6, "count": 100},
]

RNG = np.random.default_rng(12345)


def compute_baseline_counts(total_snps):
    lengths = np.array([CHROM_LENGTHS[chrom] for chrom in CHROMOSOMES], dtype=np.float64)
    proportions = lengths / lengths.sum()
    raw_counts = np.floor(proportions * total_snps).astype(int)

    remainder = total_snps - raw_counts.sum()
    if remainder > 0:
        fractional = (proportions * total_snps) - raw_counts
        order = np.argsort(fractional)[::-1]
        for idx in order[:remainder]:
            raw_counts[idx] += 1

    baseline_counts = {chrom: max(1, int(count)) for chrom, count in zip(CHROMOSOMES, raw_counts)}
    return baseline_counts


def generate_chromosome_rows(chromosome, baseline_count):
    rows = []
    chrom_length = CHROM_LENGTHS[chromosome]

    # Baseline SNPs: uniform sampling across the chromosome with uniform p-values
    positions = np.linspace(1, chrom_length, baseline_count, dtype=int)
    trait_a_p = np.clip(RNG.random(baseline_count), 1e-300, 1.0)
    trait_b_p = np.clip(RNG.random(baseline_count), 1e-300, 1.0)
    trait_a_pp = -np.log10(trait_a_p)
    trait_b_pp = -np.log10(trait_b_p)

    for idx, position in enumerate(positions):
        unique_eps = (idx + 1) * 1e-6
        rows.append(
            {
                "chr": chromosome,
                "pos": int(position),
                "values": {
                    "traitA_true-log10p": round(trait_a_pp[idx] + unique_eps, 6),
                    "traitB_true-log10p": round(trait_b_pp[idx] + 2 * unique_eps, 6),
                },
            }
        )

    # Peak SNPs mimicking LD decay
    for peak in [p for p in PEAKS if p["chrom"] == chromosome]:
        peak_center = int(peak["center_frac"] * chrom_length)
        half_window = int(peak["decay"] * 12)
        distances = np.linspace(-half_window, half_window, peak["count"])
        for idx, dist in enumerate(distances):
            position = peak_center + int(dist)
            if position < 1 or position > chrom_length:
                continue
            decay_factor = np.exp(-abs(dist) / peak["decay"])
            peak_height = peak["height"] * decay_factor
            ld_noise = RNG.normal(0, 0.02)
            unique_eps = (idx + 1) * 5e-6
            trait_a = round(max(peak_height + ld_noise, 0.0001) + unique_eps, 6)
            trait_b = round(max(peak_height * 0.96 + ld_noise, 0.0001) + 2 * unique_eps, 6)
            rows.append(
                {
                    "chr": chromosome,
                    "pos": position,
                    "values": {
                        "traitA_true-log10p": trait_a,
                        "traitB_true-log10p": trait_b,
                    },
                }
            )

    rows.sort(key=lambda record: record["pos"])
    return rows


def write_bgenie_file(path, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(BGENIE_HEADER)
        for row in rows:
            handle.write(
                f"{row['chr']} {row['pos']} {row['values']['traitA_true-log10p']} "
                f"{row['values']['traitB_true-log10p']}\n"
            )


def canonical_sort(df):
    df = df.copy()
    df["pp"] = df["pp"].round(6)
    df["original_index"] = df["original_index"].astype(int)
    df["chr"] = df["chr"].astype(int)
    df["pos"] = df["pos"].astype(int)
    return df.sort_values(
        by=["chr", "pos", "pp", "original_index"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)


def run_python_pipeline(chromosomes, tmpdir, threshold):
    shared_dir = Path(tmpdir) / "python_chunks"
    shared_dir.mkdir()
    output_dir = Path(tmpdir) / "python_output"
    output_dir.mkdir()

    sort_script = REPO_ROOT / "divide_and_conquer" / "sort.py"
    merge_script = REPO_ROOT / "divide_and_conquer" / "pmerge_sort.py"

    combined_path = Path(tmpdir) / "python_all_chromosomes.txt"
    write_combined_bgenie(combined_path, chromosomes)
    total_rows = sum(len(rows) for rows in chromosomes.values())

    subprocess.run(
        [
            sys.executable,
            str(sort_script),
            str(combined_path),
            "--shared_dir",
            str(shared_dir),
            "--pp_threshold",
            str(int(threshold)),
            "--chunksize",
            str(total_rows),
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    subprocess.run(
        [
            sys.executable,
            str(merge_script),
            "--input_dir",
            str(shared_dir),
            "--output_dir",
            str(output_dir),
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    final_csv = output_dir / "final_sorted_data.csv"
    df = pd.read_csv(final_csv)
    df = canonical_sort(df)
    df.to_csv(final_csv, index=False)
    return df, final_csv


def write_combined_bgenie(path, chromosomes):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(BGENIE_HEADER)
        for chrom in sorted(chromosomes):
            for row in chromosomes[chrom]:
                handle.write(
                    f"{row['chr']} {row['pos']} {row['values']['traitA_true-log10p']} "
                    f"{row['values']['traitB_true-log10p']}\n"
                )


def run_cpp_pipeline(chromosomes, tmpdir, threshold, chunk_rows):
    logsort_binary, logthin_binary = ensure_cpp_binaries()
    combined_path = Path(tmpdir) / "all_chromosomes.txt"
    write_combined_bgenie(combined_path, chromosomes)

    cpp_output = Path(tmpdir) / "cpp_sorted.csv"
    cpp_chunks = Path(tmpdir) / "cpp_chunks"
    cpp_chunks.mkdir()

    subprocess.run(
        [
            str(logsort_binary),
            "--input",
            str(combined_path),
            "--output",
            str(cpp_output),
            "--threshold",
            str(threshold),
            "--chunk-rows",
            str(chunk_rows),
            "--tmpdir",
            str(cpp_chunks),
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    df = pd.read_csv(cpp_output)
    df = canonical_sort(df)
    df.to_csv(cpp_output, index=False)
    return df, logthin_binary, cpp_output


def apply_python_thinning(sorted_df, factor):
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


def run_cpp_thinning(logthin_binary, sorted_csv, tmpdir, factor):
    thinned_output = Path(tmpdir) / "cpp_thinned.csv"
    subprocess.run(
        [
            str(logthin_binary),
            "--input",
            str(sorted_csv),
            "--output",
            str(thinned_output),
            "--factor",
            str(factor),
        ],
        cwd=REPO_ROOT,
        check=True,
    )
    df = pd.read_csv(thinned_output)
    df = canonicalize_thinned(df)
    df.to_csv(thinned_output, index=False)
    return df


def run_python_thinning_script(sorted_csv, factor):
    thinning_script = REPO_ROOT / "thin_sorted_pvalues.py"
    subprocess.run(
        [
            sys.executable,
            str(thinning_script),
            str(sorted_csv),
            str(factor),
            "python",
        ],
        cwd=REPO_ROOT,
        check=True,
    )
def canonicalize_thinned(df):
    df = df.copy()
    df["pp"] = df["pp"].round(6)
    df["row_number"] = df["row_number"].astype(int)
    df["original_index"] = df["original_index"].astype(int)
    df["chr"] = df["chr"].astype(int)
    df["pos"] = df["pos"].astype(int)
    return df.sort_values("row_number").reset_index(drop=True)


def assert_sorted_tables_equal(python_sorted, cpp_sorted):
    if len(python_sorted) != len(cpp_sorted):
        raise AssertionError("Sorted tables differ in length.")

    np.testing.assert_array_equal(python_sorted["chr"].to_numpy(), cpp_sorted["chr"].to_numpy())
    np.testing.assert_array_equal(python_sorted["pos"].to_numpy(), cpp_sorted["pos"].to_numpy())
    np.testing.assert_array_equal(python_sorted["original_index"].to_numpy(), cpp_sorted["original_index"].to_numpy())
    np.testing.assert_allclose(
        python_sorted["pp"].to_numpy(),
        cpp_sorted["pp"].to_numpy(),
        rtol=1e-4,
        atol=1e-5,
    )


def assert_thinned_tables_equal(python_thinned, cpp_thinned):
    merged = python_thinned.merge(
        cpp_thinned,
        on="row_number",
        how="outer",
        suffixes=("_py", "_cpp"),
        indicator=True,
    )

    if not merged["_merge"].eq("both").all():
        diff = merged[merged["_merge"] != "both"]
        raise AssertionError(f"Thinned tables mismatch: {diff.head()}")

    for column in ("chr", "pos", "original_index"):
        if not (merged[f"{column}_py"] == merged[f"{column}_cpp"]).all():
            diff = merged[merged[f"{column}_py"] != merged[f"{column}_cpp"]]
            raise AssertionError(f"Column {column} differs in thinned tables: {diff.head()}")

    np.testing.assert_allclose(
        merged["pp_py"].to_numpy(),
        merged["pp_cpp"].to_numpy(),
        rtol=1e-4,
        atol=1e-5,
    )


def test_large_pipeline_consistency():
    baseline_counts = compute_baseline_counts(BACKGROUND_SNPS)
    chromosomes = {
        chrom: generate_chromosome_rows(chrom, baseline_counts[chrom])
        for chrom in CHROMOSOMES
    }
    threshold = 0.0
    chunk_rows = 512
    thinning_factor = THINNING_FACTOR

    with tempfile.TemporaryDirectory() as tmpdir:
        python_sorted, python_sorted_path = run_python_pipeline(chromosomes, tmpdir, threshold)
        python_sorted = python_sorted.reset_index(drop=True)
        cpp_sorted, logthin_binary, cpp_sorted_path = run_cpp_pipeline(chromosomes, tmpdir, threshold, chunk_rows)
        cpp_sorted = cpp_sorted.reset_index(drop=True)

        assert_sorted_tables_equal(python_sorted, cpp_sorted)

        run_python_thinning_script(python_sorted_path, thinning_factor)
        python_thinned = apply_python_thinning(python_sorted, thinning_factor)
        cpp_thinned = run_cpp_thinning(logthin_binary, cpp_sorted_path, tmpdir, thinning_factor)
        assert_thinned_tables_equal(python_thinned, cpp_thinned)


if __name__ == "__main__":
    test_large_pipeline_consistency()
    print("Large pipeline consistency test passed.")
