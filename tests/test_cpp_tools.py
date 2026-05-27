import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import subprocess

from tests.cpp_utils import REPO_ROOT, ensure_cpp_binaries

BGENIE_HEADER = "chr pos traitA_true-log10p traitB_true-log10p\n"


def write_bgenie_file(path, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(BGENIE_HEADER)
        for row in rows:
            handle.write(
                f"{row['chr']} {row['pos']} {row['values']['traitA_true-log10p']} "
                f"{row['values']['traitB_true-log10p']}\n"
            )


def build_expected_records(chromosomes):
    records = []
    for chrom, rows in chromosomes.items():
        for column in ("traitA_true-log10p", "traitB_true-log10p"):
            for row in rows:
                records.append(
                    {
                        "chr": chrom,
                        "pos": row["pos"],
                        "pp": row["values"][column],
                    }
                )
    records.sort(key=lambda r: r["pp"], reverse=True)
    return records


def test_cpp_sort_and_thin_pipeline():
    chromosomes = {
        1: [
            {"chr": 1, "pos": 111, "values": {"traitA_true-log10p": 12.1, "traitB_true-log10p": 7.2}},
            {"chr": 1, "pos": 112, "values": {"traitA_true-log10p": 11.0, "traitB_true-log10p": 6.8}},
            {"chr": 1, "pos": 113, "values": {"traitA_true-log10p": 9.5, "traitB_true-log10p": 5.0}},
        ],
        2: [
            {"chr": 2, "pos": 211, "values": {"traitA_true-log10p": 13.5, "traitB_true-log10p": 8.1}},
            {"chr": 2, "pos": 212, "values": {"traitA_true-log10p": 10.3, "traitB_true-log10p": 6.2}},
            {"chr": 2, "pos": 213, "values": {"traitA_true-log10p": 8.7, "traitB_true-log10p": 4.4}},
        ],
    }

    expected_records = build_expected_records(chromosomes)

    logsort_binary, logthin_binary = ensure_cpp_binaries()

    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / "all_chromosomes.txt"
        write_bgenie_file(
            input_path,
            chromosomes[1] + chromosomes[2],
        )

        sorted_csv = Path(tmpdir) / "final_sorted_data.csv"
        temp_chunks = Path(tmpdir) / "chunks"
        temp_chunks.mkdir()

        subprocess.run(
            [
                str(logsort_binary),
                "--input",
                str(input_path),
                "--output",
                str(sorted_csv),
                "--threshold",
                "0",
                "--chunk-rows",
                "2",
                "--tmpdir",
                str(temp_chunks),
            ],
            cwd=REPO_ROOT,
            check=True,
        )

        assert sorted_csv.exists(), "logsort did not create the expected output CSV."

        final_df = pd.read_csv(sorted_csv)
        actual_records = final_df[["chr", "pos", "pp"]].to_dict("records")
        assert actual_records == expected_records

        thinned_csv = Path(tmpdir) / "thinned_sorted_data.csv"
        subprocess.run(
            [
                str(logthin_binary),
                "--input",
                str(sorted_csv),
                "--output",
                str(thinned_csv),
                "--factor",
                "1.5",
            ],
            cwd=REPO_ROOT,
            check=True,
        )

        assert thinned_csv.exists(), "logthin did not create the expected output CSV."

        thinned_df = pd.read_csv(thinned_csv)
        assert "row_number" in thinned_df.columns
        thinned_records = thinned_df[["chr", "pos", "pp"]].to_dict("records")
        for record in thinned_records:
            assert record in expected_records


if __name__ == "__main__":
    test_cpp_sort_and_thin_pipeline()
    print("C++ pipeline test passed.")
