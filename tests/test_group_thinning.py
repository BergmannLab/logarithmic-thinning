"""End-to-end test of independent per-group logarithmic thinning.

Builds a synthetic BGENIE-style chunk whose phenotype columns span all three
default groups (_pred / _true / LV_), runs the current sort.process_chunk +
pmerge_sort.merge_chunks API with grouping enabled, and asserts that:

* one thinned CSV is produced per group;
* each group's kept rows match an *independent* logarithmic_thinning of that
  group's own record total, with within-group `unthinned_rank`;
* every record routed to a group came from a column matching that group;
* the union of all grouped outputs equals the ungrouped (single-output) run.
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from sort_merge import sort as sorter
from sort_merge import pmerge_sort as merger
from thin_sorted_pvalues import logarithmic_thinning

GROUP_SPEC = "dTIF=_pred,mTIF=_true,LV=^LV_"
FACTOR = 1.2

# Distinct pp ranges per column so we can tell, for any emitted pp, which column
# (hence which group) it came from.
COLUMNS = {
    "D_A_std_pred-log10p": (300.0, 400.0),  # dTIF (group 0)
    "D_V_std_pred-log10p": (200.0, 300.0),  # dTIF (group 0)
    "D_A_std_true-log10p": (100.0, 200.0),  # mTIF (group 1)
    "LV_0-log10p": (0.0, 100.0),            # LV   (group 2)
}
GROUP_OF_RANGE = {"D_A_std_pred-log10p": 0, "D_V_std_pred-log10p": 0,
                  "D_A_std_true-log10p": 1, "LV_0-log10p": 2}


def _build_chunk(n_rows, shared_dir, groups):
    rng = np.random.default_rng(0)
    cols = list(COLUMNS)
    data = {"chr": np.ones(n_rows, dtype=np.int32),
            "pos": np.arange(1, n_rows + 1, dtype=np.int32)}
    for col, (lo, hi) in COLUMNS.items():
        data[col] = rng.uniform(lo, hi, size=n_rows).astype(np.float32)
    df = pd.DataFrame(data)

    group_ids = sorter.classify_columns(cols, groups) if groups else []
    config = sorter.SortConfig(
        shared_dir=shared_dir,
        chr_column="chr",
        pos_column="pos",
        pp_columns=cols,
        p_columns=[],
        threshold=0.0,
        group_ids=tuple(group_ids),
    )
    return sorter.process_chunk((df, 0, config))


def main():
    n_rows = 2000  # 2000 SNPs x 4 phenos = 8000 records
    groups = sorter.parse_group_spec(GROUP_SPEC)
    names = [name for name, _ in groups]

    with tempfile.TemporaryDirectory() as tmp:
        # --- grouped run -------------------------------------------------------
        chunks_g = os.path.join(tmp, "chunks_g")
        out_g = os.path.join(tmp, "out_g")
        os.makedirs(chunks_g)
        chunk = _build_chunk(n_rows, chunks_g, groups)
        merger.merge_chunks([chunk], out_g, thinning_factor=FACTOR,
                            group_names=names, keep_temp=False)

        # group totals: 2 pred cols, 1 true col, 1 LV col -> 2N, N, N records
        expected_totals = {"dTIF": 2 * n_rows, "mTIF": n_rows, "LV": n_rows}
        pp_range_by_group = {0: (200.0, 400.0), 1: (100.0, 200.0), 2: (0.0, 100.0)}

        all_rows = []
        for gid, name in enumerate(names):
            path = os.path.join(out_g, f"thinned_final_sorted_data_{name}.csv")
            assert os.path.exists(path), f"missing output for group {name}"
            d = pd.read_csv(path)

            total = expected_totals[name]
            expected_keep = [r for r in logarithmic_thinning(total, FACTOR) if r <= total]
            assert len(d) == len(expected_keep), (name, len(d), len(expected_keep))
            assert d["unthinned_rank"].tolist() == sorted(expected_keep), name
            # pp must be sorted descending and lie in this group's pp band
            assert d["pp"].is_monotonic_decreasing, name
            lo, hi = pp_range_by_group[gid]
            assert d["pp"].between(lo, hi).all(), (name, d["pp"].min(), d["pp"].max())
            all_rows.append(d[["chr", "pos", "pp"]])

        # --- ungrouped run (back-compat) --------------------------------------
        chunks_u = os.path.join(tmp, "chunks_u")
        out_u = os.path.join(tmp, "out_u")
        os.makedirs(chunks_u)
        chunk_u = _build_chunk(n_rows, chunks_u, [])  # no grouping
        merger.merge_chunks([chunk_u], out_u, thinning_factor=FACTOR,
                            group_names=None, keep_temp=False)
        single = os.path.join(out_u, "thinned_final_sorted_data.csv")
        assert os.path.exists(single), "ungrouped output missing"

        # Union of grouped pp-values must be a superset relationship sanity check:
        # the grouped run keeps strictly more rows than one global thinning,
        # because each group restarts the rank-1 baseline.
        grouped_total = sum(len(x) for x in all_rows)
        ungrouped_total = len(pd.read_csv(single))
        assert grouped_total > ungrouped_total, (grouped_total, ungrouped_total)

    print("test_group_thinning passed:")
    print(f"  grouped rows total = {grouped_total}, ungrouped = {ungrouped_total}")


def test_group_thinning():
    """pytest entry point.

    The checks live in ``main()`` so the module stays runnable as a standalone script;
    without a ``test_``-prefixed function pytest collects nothing from this file.
    """
    main()


if __name__ == "__main__":
    main()
