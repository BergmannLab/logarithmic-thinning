#!/usr/bin/env bash
# Full logarithmic-thinning pipeline from raw BGENIE per-chromosome outputs.
#
#   ./run_pipeline.sh <bgenie_dir> <out_dir> [thinning_factor] [glob]
#
# <bgenie_dir>      directory of raw BGENIE files (one per chromosome)
# <out_dir>         where chunks + final/thinned CSVs are written
# [thinning_factor] default 1.0003
# [glob]            filename pattern within <bgenie_dir>, default '*' (e.g. 'chr*.txt')
#
# Assumes BGENIE defaults: space-delimited, with `chr`/`pos` columns and one or
# more `*true-log10p` columns (auto-detected). Override in the sort step if not.
set -euo pipefail

IN="${1:?bgenie input dir}"
OUT="${2:?output dir}"
FACTOR="${3:-1.0003}"
GLOB="${4:-*}"
HERE="$(cd "$(dirname "$0")" && pwd)"
CHUNKS="$OUT/chunks"
mkdir -p "$CHUNKS"

# 1) sort each chromosome into a shared chunk dir (chunks accumulate across files)
for f in "$IN"/$GLOB; do
  echo "[sort] $f"
  python "$HERE/divide_and_conquer/sort.py" "$f" \
    --shared_dir "$CHUNKS" --auto-detect --pp-threshold 0
done

# 2) merge all chunks -> $OUT/final_sorted_data.csv (+ .npy)
echo "[merge]"
python "$HERE/divide_and_conquer/pmerge_sort.py" --input_dir "$CHUNKS" --output_dir "$OUT"

# 3) logarithmic thinning -> $OUT/thinned_final_sorted_data.csv
echo "[thin]"
python "$HERE/thin_sorted_pvalues.py" "$OUT/final_sorted_data.csv" \
  --thinning-factor "$FACTOR" --method python \
  --output "$OUT/thinned_final_sorted_data.csv"

echo "Done -> $OUT/thinned_final_sorted_data.csv"
