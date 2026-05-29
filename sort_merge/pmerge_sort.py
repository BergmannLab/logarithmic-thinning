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

# Packed (unaligned) record: pp f32 + indices i32 + chr i16 + pos i32 = 14 bytes.
RECORD_DTYPE = np.dtype(
    {
        "names": ["pp", "indices", "chr", "pos"],
        "formats": [np.float32, np.int32, np.int16, np.int32],
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


def _convert_to_memmap(npz_path: str, tmp_dir: str) -> str:
    """Load a compressed sorted chunk and rewrite it as a memmappable .npy.

    Returns the path to the structured ``.npy`` (sorted descending by pp, same
    order as the input). Peak memory is one chunk.
    """
    data = np.load(npz_path)
    pp = data["pp"]
    n = pp.shape[0]

    rec = np.empty(n, dtype=RECORD_DTYPE)
    rec["pp"] = pp
    rec["indices"] = data["indices"]
    rec["chr"] = data["chr"]
    rec["pos"] = data["pos"]
    del data, pp

    base = os.path.splitext(os.path.basename(npz_path))[0]
    out_path = os.path.join(tmp_dir, f"{base}.mm.npy")
    np.save(out_path, rec)
    del rec
    return out_path


def _emit_selected(
    batch: np.ndarray,
    offset: int,
    keep_ranks: np.ndarray,
    csv_handle,
) -> None:
    """Sort a gathered batch descending by pp, then append only the rows whose
    global 1-based rank falls in ``keep_ranks``.

    The records in ``batch`` occupy global ranks ``(offset, offset + len(batch)]``
    once sorted. ``keep_ranks`` is the ascending array of ranks to retain (the
    logarithmic-thinning index set). Because each retained rank uniquely
    identifies a sorted position, we never materialise the full sorted batch as
    a DataFrame -- only the (typically tiny) selected subset is formatted and
    written. This is what turns the merge from a ~770 GB CSV write into a few MB.
    """
    n = batch.shape[0]
    lo = int(np.searchsorted(keep_ranks, offset + 1, side="left"))
    hi = int(np.searchsorted(keep_ranks, offset + n, side="right"))
    if hi <= lo:
        return  # nothing in this block survives thinning

    order = np.argsort(-batch["pp"], kind="stable")
    sel_ranks = keep_ranks[lo:hi]
    # sorted position of rank r is (r - offset - 1); map back to the original
    # row in batch via the argsort order.
    sel = batch[order[sel_ranks - offset - 1]]
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


def _stream_kway_merge(
    memmaps: List[np.ndarray], csv_path: str, block: int, keep_ranks: np.ndarray
) -> int:
    """Single-pass k-way merge of sorted-descending memmaps, thinned inline.

    Each ``memmaps[k]`` is a structured array sorted descending by ``pp``. The
    merge streams the globally-sorted order but only writes rows whose global
    rank is in ``keep_ranks``; it returns the number of rows actually written.

    Invariant per iteration: we only emit records whose ``pp`` is >= ``t``,
    where ``t`` is the largest "block tail" among chunks that still have
    unloaded data. Every unloaded record is strictly < ``t``, so the gathered
    ``pp >= t`` records are guaranteed to be the next ones globally and can be
    sorted among themselves and flushed. The chunk owning ``t`` always has its
    whole loaded block emitted, guaranteeing forward progress.
    """
    k = len(memmaps)
    lengths = [m.shape[0] for m in memmaps]
    pos = [0] * k

    offset = 0  # global rank already passed (number of records consumed)
    written = 0
    with open(csv_path, "w", buffering=1 << 20, newline="") as fh:
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
            before = offset + batch.shape[0]
            lo = int(np.searchsorted(keep_ranks, offset + 1, side="left"))
            hi = int(np.searchsorted(keep_ranks, before, side="right"))
            _emit_selected(batch, offset, keep_ranks, fh)
            written += hi - lo
            offset = before

    return written


def merge_chunks(
    input_files,
    output_dir,
    thinning_factor=1.0003,
    n_workers=None,  # accepted for CLI/back-compat; merge is single-process
    tmp_dir=None,
    mem_budget_gb=50.0,
    keep_temp=False,
):
    """Merge sorted ``.npz`` chunks, thin inline, write thinned CSV.

    The merge streams the globally-sorted order but only writes the rows kept by
    logarithmic thinning (factor ``thinning_factor``), producing
    ``output_dir/thinned_final_sorted_data.csv`` directly. The full
    fully-expanded ordering (~10^10 rows) is never written to disk, since the
    downstream thin step would discard all but a few tens of thousands of rows.
    """
    input_files = list(input_files)
    if not input_files:
        raise ValueError("No input chunk files to merge.")

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
    for idx, f in enumerate(input_files, 1):
        mm_paths.append(_convert_to_memmap(f, tmp_dir))
        if idx % 10 == 0 or idx == k:
            print(f"  converted {idx}/{k} chunks to memmap", flush=True)
    print(f"Conversion done in {time.time() - t0:.1f}s.", flush=True)

    memmaps = [np.load(p, mmap_mode="r") for p in mm_paths]

    # The thinning rank-set depends only on the total record count, which is the
    # sum of chunk lengths -- known before we merge a single record.
    total_rows = int(sum(m.shape[0] for m in memmaps))
    keep_ranks = np.asarray(
        logarithmic_thinning(total_rows, thinning_factor), dtype=np.int64
    )

    csv_path = os.path.join(output_dir, "thinned_final_sorted_data.csv")
    t1 = time.time()
    total = _stream_kway_merge(memmaps, csv_path, block, keep_ranks)
    print(
        f"Wrote {total:,} thinned records to {csv_path} in {time.time() - t1:.1f}s.",
        flush=True,
    )

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

    return csv_path


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
    args = parser.parse_args()

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
    )


if __name__ == "__main__":
    main()
