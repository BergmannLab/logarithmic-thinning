#!/usr/bin/env python3
"""Check the thinning SELECTION behind a published figure -- not the full pipeline.

Figure 3a of the retinal-vasculometry GWAS paper is drawn from three thinned SNP tables
produced by this package. Each row of those tables carries the rank it held in the full
sorted p-value list, so the published selection can be compared directly against a freshly
computed one.

This compares WHICH ROWS were kept, in milliseconds. It does not regenerate the row
contents: that needs the end-to-end pipeline (sort, merge, then thin), which takes minutes
for the TIF sets and up to hours for the 1024 latent variables. Use ``run_pipeline.sh`` for
that, and see the scope note in ``tests/test_published_figure_inputs.py``.

    python examples/reproduce_figure3a_ranks.py <dir-with-*_thinned_snps.csv>

The script recovers the pre-thinning row count from the largest recorded rank, regenerates
the kept ranks with ``logarithmic_thinning``, and reports whether the two agree exactly.

Note on the file format: these CSVs were written by an early version of the pipeline that
emitted the rank column without a header name, so the first data row's rank (``1``) was
consumed as the column heading -- which is why the header reads ``...,original_index,1`` and
the first listed rank is 2. This script accounts for that.

No source data is needed to verify the algorithm itself: the kept-rank set depends only on
the row count and the thinning factor. ``tests/test_published_figure_inputs.py`` pins the
same three cases by checksum and runs anywhere, without the CSVs.
"""

import csv
import sys
from pathlib import Path

from logarithmic_thinning import logarithmic_thinning

FACTOR = 1.0003  # the pipeline default, and what produced the published tables
TRAIT_SETS = ("mtifs", "dtifs", "lvs")
RANK_COLUMN = 4  # chr, pos, pp, original_index, <rank>


def published_ranks(path):
    """Return the ranks recorded in a published thinned CSV, header quirk included."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        ranks = [int(header[RANK_COLUMN])]  # rank 1, stranded in the header
        ranks.extend(int(row[RANK_COLUMN]) for row in reader)
    return ranks


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2

    directory = Path(sys.argv[1])
    failures = 0

    for name in TRAIT_SETS:
        path = directory / f"{name}_thinned_snps.csv"
        if not path.exists():
            print(f"{name:6s} SKIP  not found: {path}")
            continue

        expected = published_ranks(path)
        total_rows = max(expected)
        regenerated = [r for r in logarithmic_thinning(total_rows, FACTOR) if 1 <= r <= total_rows]

        ok = regenerated == expected
        failures += not ok
        print(
            f"{name:6s} {'MATCH' if ok else 'DIFFER'}  "
            f"rows before thinning = {total_rows:>12,}  kept = {len(expected):>7,}"
        )
        if not ok:
            print(f"       regenerated {len(regenerated):,} ranks; first difference at "
                  f"{next(i for i, (a, b) in enumerate(zip(regenerated, expected)) if a != b)}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
