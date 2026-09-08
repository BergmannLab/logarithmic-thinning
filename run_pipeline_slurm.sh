#!/usr/bin/env bash
# SLURM variant of run_pipeline.sh.
#
# Submits a job array for the per-chromosome sort step, then a single merge
# job that depends on `afterok` of the sort array. The merge applies
# logarithmic thinning inline, so there is no separate thin job; if any sort
# task fails, the merge job is auto-cancelled.
#
# Usage:
#   ./run_pipeline_slurm.sh <bgenie_dir> <out_dir> \
#                           [thinning_factor=1.0003] \
#                           [glob='*'] \
#                           [chunksize=500000] \
#                           [col_pattern='-log10p$'] \
#                           [groups='dTIF=_pred,mTIF=_true,LV=^LV_']
#
# Phenotype columns are partitioned into independently-thinned groups by the
# <groups> spec (comma-separated 'name=regex' pairs; each column joins the first
# group whose regex matches, a column matching none aborts the sort task). The
# default groups dTIFs (_pred), mTIFs (_true) and LVs (^LV_); pass '' to thin
# all columns together.
#
# Layout written under <out_dir>:
#   chunks/                   -- intermediate .npz chunks (one per sort task)
#   slurm/inputs.txt          -- one input file path per line (array index)
#   slurm/pp_columns.txt      -- the comma-separated --pp-columns list
#   slurm/sort.sbatch         -- generated array job script
#   slurm/merge.sbatch        -- generated merge (+ inline thin) job script
#   slurm/logs/               -- per-job stdout/stderr
#   thinned_final_sorted_data_<group>.csv  -- one thinned output per group
#                                             (thinned_final_sorted_data.csv when
#                                             ungrouped)
#
# Resources default to:
#   sort  : 8 CPUs, 200G, 2h  (account/partition from the environment,
#                            see the site-configuration block below)
#   merge : 1 CPU,  96G,  4h
# Sort memory is generous because each worker pickles its full chunk
# DataFrame (~5 GB for 500k rows × 1000 phenos × float32) and we may have
# 3+ workers running in parallel on the largest chromosomes.
# Merge is a single-process streaming k-way merge whose peak RAM is bounded by
# --mem-budget-gb (set below), independent of total dataset size; 96G leaves
# ample headroom over the 48G budget for memmap page cache and pandas overhead.
# Thinning is fused into the merge: the globally-sorted stream is produced in
# rank order, and only the logarithmically-spaced kept ranks are ever written,
# so the merge emits ~tens of thousands of rows instead of a ~770 GB CSV.
# Edit the SBATCH_* variables below to change them.
set -euo pipefail

IN="${1:?bgenie input dir}"
OUT="${2:?output dir}"
FACTOR="${3:-1.0003}"
GLOB="${4:-*}"
CHUNK="${5:-500000}"
PATTERN="${6:-}"
[ -z "$PATTERN" ] && PATTERN='-log10p$'
# 7th arg may be '' to disable grouping; only fall back to the default when the
# arg is entirely absent rather than an explicit empty string.
GROUP_SPEC="${7-dTIF=_pred,mTIF=_true,LV=^LV_}"

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$(cd "$(dirname "$OUT")" && pwd)/$(basename "$OUT")"  # absolutise
CHUNKS="$OUT/chunks"
SLURM_DIR="$OUT/slurm"
LOG_DIR="$SLURM_DIR/logs"

# ---------------------------------------------------------------------------
# Site configuration. Every value is empty by default and overridden from the
# environment, so the driver runs on any SLURM cluster out of the box:
#
#   SBATCH_ACCOUNT=my_account SBATCH_PARTITION=my_partition ./run_pipeline_slurm.sh ...
#
# An empty value means "omit the directive" rather than "pass an empty one" --
# clusters that need no account, or that have a default partition, need nothing
# set. SBATCH_MODULES is any shell prologue the job must run before `python` is
# on PATH (module loads, a venv activation, `pixi shell-hook`, ...).
# ---------------------------------------------------------------------------
SBATCH_ACCOUNT="${SBATCH_ACCOUNT:-}"
SBATCH_PARTITION="${SBATCH_PARTITION:-}"
SBATCH_MODULES="${SBATCH_MODULES:-}"

# Assembled once; interpolated into both generated sbatch scripts. Each entry
# carries its own trailing newline so an unset value leaves no blank line
# between the shebang and the remaining directives.
# `if` rather than `[ ... ] && ...`: this runs under `set -e`, where a bare test
# that fails is the script's exit status and would abort the run whenever a value
# is legitimately unset.
SBATCH_SITE=""
if [ -n "$SBATCH_ACCOUNT" ]; then
  SBATCH_SITE="${SBATCH_SITE}#SBATCH --account=${SBATCH_ACCOUNT}
"
fi
if [ -n "$SBATCH_PARTITION" ]; then
  SBATCH_SITE="${SBATCH_SITE}#SBATCH --partition=${SBATCH_PARTITION}
"
fi

[ -d "$IN" ] || { echo "ERROR: input dir does not exist: $IN" >&2; exit 1; }
if [ -d "$CHUNKS" ] && [ -n "$(ls -A "$CHUNKS" 2>/dev/null)" ]; then
  echo "ERROR: $CHUNKS already exists and is non-empty." >&2
  echo "Remove it (rm -rf '$CHUNKS') or choose a different <out_dir>." >&2
  exit 1
fi
command -v sbatch >/dev/null || { echo "ERROR: sbatch not on PATH." >&2; exit 1; }
mkdir -p "$CHUNKS" "$SLURM_DIR" "$LOG_DIR"

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

# inputs.txt: one absolute path per line; line N == array task N.
INPUTS="$SLURM_DIR/inputs.txt"
: > "$INPUTS"
for f in "${files[@]}"; do
  printf '%s\n' "$f" >> "$INPUTS"
done
N="${#files[@]}"

# Derive --pp-columns from the first file's header.
PP_COLS=$(read_header "${files[0]}" | tr ' \t' '\n\n' | grep -E -e "$PATTERN" | paste -sd,)
if [ -z "$PP_COLS" ]; then
  echo "ERROR: no header column matches /$PATTERN/ in ${files[0]}" >&2; exit 1
fi
N_COLS=$(echo "$PP_COLS" | tr ',' '\n' | wc -l)

# Sanity-check: compare against the last file's matching columns.
if [ "$N" -gt 1 ]; then
  LAST_COLS=$(read_header "${files[$N-1]}" | tr ' \t' '\n\n' | grep -E -e "$PATTERN" | paste -sd,)
  if [ "$PP_COLS" != "$LAST_COLS" ]; then
    echo "ERROR: header columns differ between ${files[0]} and ${files[$N-1]}." >&2
    echo "       Aborting before any jobs are submitted." >&2
    exit 1
  fi
fi

COLS_FILE="$SLURM_DIR/pp_columns.txt"
printf '%s\n' "$PP_COLS" > "$COLS_FILE"

echo "Submitting SLURM pipeline:"
echo "  inputs       : $N file(s)"
echo "  phenotypes   : $N_COLS column(s) matching /$PATTERN/"
echo "  groups       : ${GROUP_SPEC:-<none, single output>}"
echo "  chunks dir   : $CHUNKS"
echo "  slurm dir    : $SLURM_DIR"

# ---------------------------------------------------------------------------
# Generate the three sbatch scripts. Substitutions like $OUT/$HERE/$CHUNK are
# expanded NOW (script generation). Anything that must be resolved at job
# run-time is escaped (\$).
# ---------------------------------------------------------------------------

SORT_SBATCH="$SLURM_DIR/sort.sbatch"
cat > "$SORT_SBATCH" <<EOF
#!/usr/bin/env bash
${SBATCH_SITE}#SBATCH --job-name=logthin-sort
#SBATCH --output=$LOG_DIR/sort-%A_%a.out
#SBATCH --error=$LOG_DIR/sort-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=02:00:00
set -euo pipefail
$SBATCH_MODULES

FILE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$INPUTS")
COLS=\$(cat "$COLS_FILE")
echo "[\$(date -Is)] sort task \$SLURM_ARRAY_TASK_ID -> \$FILE"

python "$HERE/sort_merge/sort.py" "\$FILE" \\
  --shared_dir "$CHUNKS" \\
  --pp-columns "\$COLS" \\
  --pp-threshold 0 \\
  --chunksize $CHUNK \\
  --groups "$GROUP_SPEC" \\
  --workers \$SLURM_CPUS_PER_TASK
EOF

MERGE_SBATCH="$SLURM_DIR/merge.sbatch"
cat > "$MERGE_SBATCH" <<EOF
#!/usr/bin/env bash
${SBATCH_SITE}#SBATCH --job-name=logthin-merge
#SBATCH --output=$LOG_DIR/merge-%j.out
#SBATCH --error=$LOG_DIR/merge-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=04:00:00
set -euo pipefail
$SBATCH_MODULES

echo "[\$(date -Is)] merge starting"
python "$HERE/sort_merge/pmerge_sort.py" \\
  --input_dir "$CHUNKS" \\
  --output_dir "$OUT" \\
  --thinning-factor $FACTOR \\
  --groups "$GROUP_SPEC" \\
  --mem-budget-gb 48
EOF

# ---------------------------------------------------------------------------
# Submit the DAG. Thinning is fused into the merge, so there is no separate
# thin job: the merge writes thinned_final_sorted_data.csv directly.
# ---------------------------------------------------------------------------
SORT_ID=$(sbatch --parsable --array=1-$N "$SORT_SBATCH")
MERGE_ID=$(sbatch --parsable --dependency=afterok:$SORT_ID "$MERGE_SBATCH")

echo
echo "Submitted:"
echo "  sort  : $SORT_ID (array 1-$N)"
echo "  merge : $MERGE_ID (afterok:$SORT_ID)  [merges + thins inline]"
echo
echo "Logs:    $LOG_DIR/"
if [ -n "$GROUP_SPEC" ]; then
  echo "Output:  $OUT/thinned_final_sorted_data_<group>.csv (one per group, once merge succeeds)"
else
  echo "Output:  $OUT/thinned_final_sorted_data.csv (once merge succeeds)"
fi
echo "Watch:   squeue -u \$USER -j $SORT_ID,$MERGE_ID"
