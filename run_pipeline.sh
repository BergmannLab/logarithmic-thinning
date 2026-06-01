#!/usr/bin/env bash
# Full logarithmic-thinning pipeline from raw BGenie per-chromosome outputs.
#
# Usage:
#   ./run_pipeline.sh <bgenie_dir> <out_dir> \
#                     [thinning_factor=1.0003] \
#                     [glob='*'] \
#                     [chunksize=500000] \
#                     [col_pattern='-log10p$'] \
#                     [groups='dTIF=_pred,mTIF=_true,LV=^LV_']
#
# Behaviour:
#   - For each input file matching <bgenie_dir>/<glob>, derives the list of
#     phenotype columns from the header (every column name matching
#     <col_pattern>, default = every '-log10p' column) and passes it to
#     sort_merge/sort.py as --pp-columns. The sort step writes
#     sorted .npz chunks into <out_dir>/chunks/.
#   - Phenotype columns are partitioned into independently-thinned groups by the
#     <groups> spec (comma-separated 'name=regex' pairs; each column joins the
#     first group whose regex matches it, a column matching none is an error).
#     The default groups dTIFs (_pred), mTIFs (_true) and LVs (^LV_). Pass an
#     empty string ('') to thin all columns together as one group.
#   - pmerge_sort.py then merges every chunk into a single globally-sorted
#     stream and applies logarithmic thinning inline *per group*, writing only
#     the kept rows to <out_dir>/thinned_final_sorted_data_<group>.csv (or
#     thinned_final_sorted_data.csv when ungrouped). The full (unthinned)
#     ordering is never materialised on disk.
#
# Notes:
#   - Inputs are assumed space-delimited with `chr`/`pos` columns and one or
#     more `-log10p` columns already on the -log10 scale.
#   - Plain-text and gzipped (.gz) inputs are both supported; the header is
#     read with zcat for .gz files, and pandas reads .gz natively.
#   - <out_dir>/chunks must not exist (or must be empty); the script refuses
#     to merge into a chunk dir that already has data, to avoid silently
#     mixing runs.
#   - Memory of the sort step scales roughly as
#         chunksize x n_phenotypes x 4 bytes
#     For ~1000 phenotypes and chunksize=500000 that is ~2 GB just for the
#     pp matrix, plus pandas overhead. Lower chunksize for wider inputs.
#   - To restrict the set of thinned columns, pass an extended-regex as the
#     6th argument, e.g.:
#         '_true-log10p$'             only measured-trait columns
#         'LV_.*-log10p$'             only the LV columns
#         '(LV_.*|.*_true)-log10p$'   LV + measured
set -euo pipefail

IN="${1:?bgenie input dir}"
OUT="${2:?output dir}"
FACTOR="${3:-1.0003}"
GLOB="${4:-*}"
CHUNK="${5:-500000}"
PATTERN="${6:-}"
[ -z "$PATTERN" ] && PATTERN='-log10p$'
# 7th arg may be '' to disable grouping; only fall back to the default when the
# arg is entirely absent ('-' guard) rather than an explicit empty string.
GROUP_SPEC="${7-dTIF=_pred,mTIF=_true,LV=^LV_}"

HERE="$(cd "$(dirname "$0")" && pwd)"
CHUNKS="$OUT/chunks"

[ -d "$IN" ] || { echo "ERROR: input dir does not exist: $IN" >&2; exit 1; }
if [ -d "$CHUNKS" ] && [ -n "$(ls -A "$CHUNKS" 2>/dev/null)" ]; then
  echo "ERROR: $CHUNKS already exists and is non-empty." >&2
  echo "Remove it (rm -rf '$CHUNKS') or choose a different <out_dir>." >&2
  exit 1
fi
mkdir -p "$CHUNKS"

read_header() {
  case "$1" in
    *.gz) zcat "$1" | head -1 ;;
    *)    head -1 "$1" ;;
  esac
}

shopt -s nullglob
files=("$IN"/$GLOB)
shopt -u nullglob
if [ "${#files[@]}" -eq 0 ]; then
  echo "ERROR: no files match $IN/$GLOB" >&2; exit 1
fi

# 1) sort each chromosome into the shared chunk dir
for f in "${files[@]}"; do
  COLS=$(read_header "$f" | tr ' \t' '\n\n' | grep -E -e "$PATTERN" | paste -sd,)
  if [ -z "$COLS" ]; then
    echo "ERROR: no header column matches /$PATTERN/ in $f" >&2; exit 1
  fi
  N=$(echo "$COLS" | tr ',' '\n' | wc -l)
  echo "[sort] $f  ($N phenotypes, chunksize=$CHUNK, groups='${GROUP_SPEC:-<none>}')"
  python "$HERE/sort_merge/sort.py" "$f" \
    --shared_dir "$CHUNKS" \
    --pp-columns "$COLS" \
    --pp-threshold 0 \
    --chunksize "$CHUNK" \
    --groups "$GROUP_SPEC"
done

# 2) merge all chunks, thinning inline per group -> $OUT/thinned_final_sorted_data*.csv
echo "[merge+thin]"
python "$HERE/sort_merge/pmerge_sort.py" \
  --input_dir "$CHUNKS" --output_dir "$OUT" \
  --thinning-factor "$FACTOR" \
  --groups "$GROUP_SPEC"

if [ -n "$GROUP_SPEC" ]; then
  echo "Done -> $OUT/thinned_final_sorted_data_<group>.csv"
else
  echo "Done -> $OUT/thinned_final_sorted_data.csv"
fi
