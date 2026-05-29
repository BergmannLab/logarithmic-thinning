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


def _emit_batch(batch: np.ndarray, csv_handle, write_header: bool) -> None:
    """Sort a gathered batch descending by pp and append it to the CSV."""
    order = np.argsort(-batch["pp"], kind="stable")
    batch = batch[order]
    df = pd.DataFrame(
        {
            "chr": batch["chr"],
            "pos": batch["pos"],
            "pp": batch["pp"],
            "original_index": batch["indices"],
        }
    )
    df.to_csv(csv_handle, header=write_header, index=False)


def _stream_kway_merge(
    memmaps: List[np.ndarray], csv_path: str, block: int
) -> int:
    """Single-pass k-way merge of sorted-descending memmaps into a CSV.

    Each ``memmaps[k]`` is a structured array sorted descending by ``pp``.
    Returns the total number of records written.

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

    total_written = 0
    write_header = True
    with open(csv_path, "w", buffering=1 << 20, newline="") as fh:
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
            _emit_batch(batch, fh, write_header)
            write_header = False
            total_written += batch.shape[0]

    return total_written


def merge_chunks(
    input_files,
    output_dir,
    n_workers=None,  # accepted for CLI/back-compat; merge is single-process
    tmp_dir=None,
    mem_budget_gb=50.0,
    keep_temp=False,
):
    """Merge sorted ``.npz`` chunks into ``output_dir/final_sorted_data.csv``."""
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

    csv_path = os.path.join(output_dir, "final_sorted_data.csv")
    t1 = time.time()
    total = _stream_kway_merge(memmaps, csv_path, block)
    print(
        f"Wrote {total:,} sorted records to {csv_path} in {time.time() - t1:.1f}s.",
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
    parser = argparse.ArgumentParser(description="Streaming k-way merge of sorted chunks.")
    parser.add_argument("--input_dir", required=True, help="Directory containing sorted chunk .npz files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save final_sorted_data.csv.")
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
        n_workers=args.workers,
        tmp_dir=args.tmp_dir,
        mem_budget_gb=args.mem_budget_gb,
        keep_temp=args.keep_temp,
    )


if __name__ == "__main__":
    main()
