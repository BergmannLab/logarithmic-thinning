#!/usr/bin/env python3
"""
Memory-bounded merge of pre-sorted GWAS chunks.

Each input ``.npz`` chunk (written by ``sort.py``) is already sorted in
descending ``pp`` (-log10 p) order. This module merges them all into a single
``final_sorted_data.csv``, sorted descending by ``pp``, using a *streaming
k-way merge* whose peak memory is bounded by ``num_chunks * block_size`` rather
than by the total dataset size.

Why not the old pairwise concat+argsort?  With ``--pp-threshold 0`` the dataset
is fully expanded (~N_SNPs x N_phenotypes entries, ~10^10 rows / ~140 GB at
14 bytes/row). A pairwise re-sort holds several copies of an input at once
(concatenated arrays + an int64 ``argsort`` order array that is *double* the pp
size + the sorted copies), and the final rounds merge ~the whole dataset, so
peak RAM scales with the data and OOM-kills regardless of the memory bump.

Strategy
--------
1. Convert each compressed ``.npz`` chunk to an uncompressed, memmappable
   structured ``.npy`` (one bounded pass, one chunk in RAM at a time). The
   ``.npz`` files are kept so a failed merge can be retried without re-sorting.
2. Open every converted chunk as a read-only memmap (no RAM cost) and run a
   single-pass k-way merge that emits the globally-sorted stream in blocks,
   appending directly to the output CSV.

The merge never materialises a full-length index array and never holds more
than a couple of blocks per chunk in RAM, so it fits in modest memory no matter
how large the dataset is.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from typing import List

import numpy as np
import pandas as pd

# Reuse the single source of truth for the thinning rank-set so the fused merge
# and the standalone thin_sorted_pvalues.py can never drift apart.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from thin_sorted_pvalues import UNTHINNED_COLUMN, logarithmic_thinning
from sort_merge.sort import parse_group_spec


def _output_name(group_name: str) -> str:
    """CSV filename for a group ('' => the single ungrouped output)."""
    if group_name:
        return f"thinned_final_sorted_data_{group_name}.csv"
    return "thinned_final_sorted_data.csv"

# Packed (unaligned) record: pp f32 + indices i32 + chr i16 + pos i32 + grp i8
# = 15 bytes. ``grp`` is the 0-based id of the independently-thinned phenotype
# group (always 0 when no grouping was requested at sort time).
RECORD_DTYPE = np.dtype(
    {
        "names": ["pp", "indices", "chr", "pos", "grp"],
        "formats": [np.float32, np.int32, np.int16, np.int32, np.int8],
    },
    align=False,
)

# Rough bytes/record consumed at peak during a merge step: gathered batch
# (RECORD_DTYPE = 14) + its int64 argsort order (8) + transient CSV/DataFrame
# formatting overhead. Used only to derive a safe block size from a memory
# budget; deliberately conservative.
_BYTES_PER_RECORD_PEAK = 64

# Never use a block larger than this many records per chunk, to keep iteration
# granularity reasonable even when few chunks are present.
_MAX_BLOCK = 50_000_000


def _convert_to_memmap(npz_path: str, tmp_dir: str, n_groups: int):
    """Load a compressed sorted chunk and rewrite it as a memmappable .npy.

    Returns ``(out_path, group_counts)`` where ``out_path`` is the structured
    ``.npy`` (sorted descending by pp, same order as the input) and
    ``group_counts`` is an int64 array of length ``n_groups`` counting records
    per group in this chunk (so the merge knows each group's total before
    streaming). Chunks written before grouping existed have no ``grp`` field and
    are treated as all group 0. Peak memory is one chunk.
    """
    data = np.load(npz_path)
    pp = data["pp"]
    n = pp.shape[0]

    rec = np.empty(n, dtype=RECORD_DTYPE)
    rec["pp"] = pp
    rec["indices"] = data["indices"]
    rec["chr"] = data["chr"]
    rec["pos"] = data["pos"]
    if "grp" in data.files:
        rec["grp"] = data["grp"]
    else:
        rec["grp"] = 0
    group_counts = np.bincount(rec["grp"], minlength=n_groups).astype(np.int64)
    del data, pp

    base = os.path.splitext(os.path.basename(npz_path))[0]
    out_path = os.path.join(tmp_dir, f"{base}.mm.npy")
    np.save(out_path, rec)
    del rec
    return out_path, group_counts


def _write_rows(csv_handle, sel: np.ndarray, sel_ranks: np.ndarray) -> None:
    """Append the selected records to one group's CSV (no header)."""
    df = pd.DataFrame(
        {
            "chr": sel["chr"],
            "pos": sel["pos"],
            "pp": sel["pp"],
            "original_index": sel["indices"],
            UNTHINNED_COLUMN: sel_ranks,
        }
    )
    df.to_csv(csv_handle, header=False, index=False)


def _emit_groups(
    batch: np.ndarray,
    group_offset: np.ndarray,
    keep_ranks_list: List[np.ndarray],
    handles: List,
    written: List[int],
) -> None:
    """Sort a gathered batch descending by pp, then append -- per group -- only
    the rows whose *within-group* 1-based rank falls in that group's keep set.

    ``group_offset[g]`` is the number of group-g records already emitted (their
    rank base) and is advanced in place. Since the global stream is descending
    in pp, each group's records form a descending subsequence, so the records of
    group ``g`` in this sorted batch occupy within-group ranks
    ``(group_offset[g], group_offset[g] + count_g]``. Each retained rank uniquely
    identifies a position, so only the (tiny) selected subset is ever formatted.
    """
    n = batch.shape[0]
    order = np.argsort(-batch["pp"], kind="stable")
    n_groups = len(handles)

    if n_groups == 1:
        # Fast path: every record is group 0; sorted positions are 0..n-1.
        kr = keep_ranks_list[0]
        base = int(group_offset[0])
        lo = int(np.searchsorted(kr, base + 1, side="left"))
        hi = int(np.searchsorted(kr, base + n, side="right"))
        if hi > lo:
            sel_ranks = kr[lo:hi]
            sel = batch[order[sel_ranks - base - 1]]
            _write_rows(handles[0], sel, sel_ranks)
            written[0] += hi - lo
        group_offset[0] += n
        return

    sorted_grp = batch["grp"][order]
    for g in range(n_groups):
        gpos = np.nonzero(sorted_grp == g)[0]  # ascending sorted positions
        cnt = gpos.size
        if cnt == 0:
            continue
        kr = keep_ranks_list[g]
        base = int(group_offset[g])
        lo = int(np.searchsorted(kr, base + 1, side="left"))
        hi = int(np.searchsorted(kr, base + cnt, side="right"))
        if hi > lo:
            sel_ranks = kr[lo:hi]
            # within-group rank r -> index (r - base - 1) into this group's
            # sorted subsequence (gpos) -> sorted-batch position -> original row.
            sel = batch[order[gpos[sel_ranks - base - 1]]]
            _write_rows(handles[g], sel, sel_ranks)
            written[g] += hi - lo
        group_offset[g] += cnt


def _stream_kway_merge(
    memmaps: List[np.ndarray],
    csv_paths: List[str],
    block: int,
    keep_ranks_list: List[np.ndarray],
) -> List[int]:
    """Single-pass k-way merge of sorted-descending memmaps, thinned inline.

    Each ``memmaps[k]`` is a structured array sorted descending by ``pp``. The
    merge streams the globally-sorted order but only writes rows whose
    within-group rank is in that group's ``keep_ranks_list`` entry, fanning out
    to one CSV per group. Returns the number of rows written per group.

    Invariant per iteration: we only emit records whose ``pp`` is >= ``t``,
    where ``t`` is the largest "block tail" among chunks that still have
    unloaded data. Every unloaded record is strictly < ``t``, so the gathered
    ``pp >= t`` records are guaranteed to be the next ones globally (hence the
    next within every group too) and can be sorted among themselves and flushed.
    The chunk owning ``t`` always has its whole loaded block emitted,
    guaranteeing forward progress.
    """
    k = len(memmaps)
    n_groups = len(csv_paths)
    lengths = [m.shape[0] for m in memmaps]
    pos = [0] * k

    group_offset = np.zeros(n_groups, dtype=np.int64)
    written = [0] * n_groups
    handles = [open(p, "w", buffering=1 << 20, newline="") for p in csv_paths]
    try:
        for fh in handles:
            fh.write("chr,pos,pp,original_index," + UNTHINNED_COLUMN + "\n")
        while True:
            active = [i for i in range(k) if pos[i] < lengths[i]]
            if not active:
                break

            # Determine the safe floor t = max tail pp among chunks that still
            # have records beyond their current block. Chunks fully covered by
            # their current block impose no floor (no unloaded data).
            t = -np.inf
            bends = {}
            for i in active:
                bend = min(pos[i] + block, lengths[i])
                bends[i] = bend
                if bend < lengths[i]:
                    tail_pp = float(memmaps[i]["pp"][bend - 1])
                    if tail_pp > t:
                        t = tail_pp

            # Gather, from each active chunk's current block, the prefix with
            # pp >= t (these are guaranteed to be the next records globally).
            parts = []
            for i in active:
                blk = memmaps[i][pos[i] : bends[i]]
                blk_pp = blk["pp"]
                if t == -np.inf:
                    count = blk_pp.shape[0]
                else:
                    # blk_pp is descending; count of values >= t equals the
                    # count of (-blk_pp) values <= -t.
                    count = int(np.searchsorted(-blk_pp, -t, side="right"))
                if count > 0:
                    parts.append(np.array(blk[:count]))  # copy out of memmap
                    pos[i] += count

            batch = parts[0] if len(parts) == 1 else np.concatenate(parts)
            _emit_groups(batch, group_offset, keep_ranks_list, handles, written)
    finally:
        for fh in handles:
            fh.close()

    return written


def merge_chunks(
    input_files,
    output_dir,
    thinning_factor=1.0003,
    n_workers=None,  # accepted for CLI/back-compat; merge is single-process
    tmp_dir=None,
    mem_budget_gb=50.0,
    keep_temp=False,
    group_names=None,
):
    """Merge sorted ``.npz`` chunks, thin inline, write thinned CSV(s).

    The merge streams the globally-sorted order but only writes the rows kept by
    logarithmic thinning (factor ``thinning_factor``). With ``group_names`` set
    (an ordered list aligned with the ``grp`` ids written by ``sort.py``), each
    group is thinned independently against its own record total and written to
    ``output_dir/thinned_final_sorted_data_<name>.csv``. Without it, all records
    form one group written to ``thinned_final_sorted_data.csv`` (legacy).
    The full fully-expanded ordering (~10^10 rows) is never written to disk.

    Returns the list of output CSV paths (one per group).
    """
    input_files = list(input_files)
    if not input_files:
        raise ValueError("No input chunk files to merge.")

    names = list(group_names) if group_names else [""]
    n_groups = len(names)

    os.makedirs(output_dir, exist_ok=True)
    if tmp_dir is None:
        tmp_dir = os.path.join(output_dir, "merge_tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    k = len(input_files)
    block = int(mem_budget_gb * (1 << 30) / max(1, k) / _BYTES_PER_RECORD_PEAK)
    block = max(100_000, min(block, _MAX_BLOCK))
    print(
        f"Merging {k} chunk(s); block={block:,} records/chunk "
        f"(mem budget {mem_budget_gb:g} GB).",
        flush=True,
    )

    t0 = time.time()
    mm_paths = []
    group_totals = np.zeros(n_groups, dtype=np.int64)
    for idx, f in enumerate(input_files, 1):
        path, counts = _convert_to_memmap(f, tmp_dir, n_groups)
        mm_paths.append(path)
        group_totals += counts
        if idx % 10 == 0 or idx == k:
            print(f"  converted {idx}/{k} chunks to memmap", flush=True)
    print(f"Conversion done in {time.time() - t0:.1f}s.", flush=True)

    memmaps = [np.load(p, mmap_mode="r") for p in mm_paths]

    # The thinning rank-set for each group depends only on that group's record
    # count -- the sum of per-chunk group counts, known before merging a record.
    keep_ranks_list: List[np.ndarray] = []
    for g, name in enumerate(names):
        total = int(group_totals[g])
        label = f" [{name}]" if name else ""
        if total == 0:
            print(f"Group{label or ' (all)'} has no records; writing header only.", flush=True)
            keep_ranks_list.append(np.empty(0, dtype=np.int64))
        else:
            print(f"Thinning group{label or ' (all)'}: {total:,} records.", flush=True)
            keep_ranks_list.append(
                np.asarray(logarithmic_thinning(total, thinning_factor), dtype=np.int64)
            )

    csv_paths = [os.path.join(output_dir, _output_name(name)) for name in names]
    t1 = time.time()
    written = _stream_kway_merge(memmaps, csv_paths, block, keep_ranks_list)
    for name, path, w in zip(names, csv_paths, written):
        label = f" [{name}]" if name else ""
        print(f"Wrote {w:,} thinned records{label} to {path}.", flush=True)
    print(f"Merge done in {time.time() - t1:.1f}s.", flush=True)

    del memmaps
    if not keep_temp:
        for p in mm_paths:
            try:
                os.remove(p)
            except OSError:
                pass
        try:
            os.rmdir(tmp_dir)
        except OSError:
            pass

    return csv_paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Streaming k-way merge of sorted chunks with inline logarithmic thinning."
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing sorted chunk .npz files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save thinned_final_sorted_data.csv.")
    parser.add_argument(
        "--thinning-factor",
        type=float,
        default=1.0003,
        help="Logarithmic thinning factor applied inline during the merge (default: 1.0003).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Accepted for back-compat; the streaming merge is single-process.",
    )
    parser.add_argument(
        "--tmp-dir",
        default=None,
        help="Scratch dir for memmapped chunks (default: <output_dir>/merge_tmp).",
    )
    parser.add_argument(
        "--mem-budget-gb",
        type=float,
        default=50.0,
        help="Approximate RAM budget used to size the per-chunk merge block (default: 50).",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the intermediate memmapped .npy files instead of deleting them.",
    )
    parser.add_argument(
        "--groups",
        default="",
        help=(
            "Comma-separated 'name=regex' pairs (the SAME spec passed to "
            "sort.py) identifying the independently-thinned phenotype groups. "
            "Only the names and their order are used here (id = position); the "
            "regexes are ignored. Each group is written to "
            "thinned_final_sorted_data_<name>.csv. Omit for a single output."
        ),
    )
    args = parser.parse_args()

    group_names = [name for name, _ in parse_group_spec(args.groups)]

    input_files = [
        os.path.join(args.input_dir, f)
        for f in os.listdir(args.input_dir)
        if f.endswith(".npz")
    ]
    if not input_files:
        print(f"ERROR: no .npz chunk files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(input_files)} chunk files to merge.", flush=True)
    merge_chunks(
        input_files,
        args.output_dir,
        thinning_factor=args.thinning_factor,
        n_workers=args.workers,
        tmp_dir=args.tmp_dir,
        mem_budget_gb=args.mem_budget_gb,
        keep_temp=args.keep_temp,
        group_names=group_names,
    )


if __name__ == "__main__":
    main()
