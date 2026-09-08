import os
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd

from sort_merge import sort as sorter
from sort_merge import pmerge_sort as merger
from thin_sorted_pvalues import logarithmic_thinning


BGENIE_HEADER = "chr pos traitA_true-log10p traitB_true-log10p\n"


def write_bgenie_file(path, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(BGENIE_HEADER)
        for row in rows:
            handle.write(
                f"{row['chr']} {row['pos']} {row['values']['traitA_true-log10p']} {row['values']['traitB_true-log10p']}\n"
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


def main():
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

    with tempfile.TemporaryDirectory() as tmpdir:
        raw_dir = os.path.join(tmpdir, "raw")
        os.makedirs(raw_dir)
        shared_dir = os.path.join(tmpdir, "chunks")
        os.makedirs(shared_dir)
        output_dir = os.path.join(tmpdir, "merged")
        os.makedirs(output_dir)

        chunk_files = []

        for chrom, rows in chromosomes.items():
            input_path = os.path.join(raw_dir, f"chr{chrom}.txt")
            write_bgenie_file(input_path, rows)

            log10p_columns = sorter.auto_detect_pp_columns(input_path, " ")
            chunk_df = pd.read_csv(
                input_path,
                sep=" ",
                usecols=["chr", "pos", *log10p_columns],
                dtype={"chr": "int8", "pos": "int32", **{col: "float32" for col in log10p_columns}},
            )

            # process_chunk takes (dataframe, start_row, SortConfig).
            config = sorter.SortConfig(
                shared_dir=shared_dir,
                chr_column="chr",
                pos_column="pos",
                pp_columns=log10p_columns,
                p_columns=(),
                threshold=0.0,
            )
            chunk_file = sorter.process_chunk((chunk_df, 0, config))
            chunk_files.append(chunk_file)

        # The merge fuses thinning in and never writes the full expanded ordering. This
        # test checks the *sorted* result, so pick a factor whose unthinned-head
        # threshold, factor/(factor - 1) ~ 2001, exceeds this tiny fixture's record
        # count -- every row is then kept.
        merger.merge_chunks(chunk_files, output_dir, thinning_factor=1.0005)

        final_csv = os.path.join(output_dir, "thinned_final_sorted_data.csv")
        final_df = pd.read_csv(final_csv)

        actual_records = final_df[["chr", "pos", "pp"]].to_dict("records")
        assert actual_records == expected_records, (
            "Merged results do not match expectation", actual_records, expected_records
        )

        thinning_factor = 1.5
        rows_to_keep = logarithmic_thinning(len(final_df), thinning_factor)
        valid_rows = [row for row in rows_to_keep if row <= len(final_df)]
        thinned_df = final_df.iloc[[row - 1 for row in valid_rows]]
        thinned_records = thinned_df[["chr", "pos", "pp"]].to_dict("records")

        for record in thinned_records:
            assert record in expected_records, ("Thinning produced an unexpected record", record)

        print("Minimal pipeline test passed. Final thinned records:")
        for record in thinned_records:
            print(record)


def test_minimal_pipeline():
    """pytest entry point.

    The checks live in ``main()`` so the module stays runnable as a standalone script;
    without a ``test_``-prefixed function pytest collects nothing from this file.
    """
    main()


if __name__ == "__main__":
    main()
