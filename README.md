# Logarithmic Thinning for the RETFound Project

Logarithmic thinning of GWAS p-values from BGenie per-chromosome outputs. Keeps
every most-significant SNP and a geometric sample of the rest, so a single
manageable file can drive Manhattan / QQ plots without losing the top of the
distribution.

## Quick start

```bash
./run_pipeline.sh <bgenie_dir> <out_dir> \
                  [thinning_factor=1.0003] \
                  [glob='*'] \
                  [chunksize=500000] \
                  [col_pattern='-log10p$'] \
                  [groups='dTIF=_pred,mTIF=_true,LV=^LV_']
```

Three stages chained:

1. **sort** — for every input file matching `<bgenie_dir>/<glob>`, derives the
   list of phenotype columns from the header (every column name matching
   `col_pattern`, default = every `-log10p` column) and passes it to
   `sort_merge/sort.py` as `--pp-columns`. Sorted `.npz` chunks land in
   `<out_dir>/chunks/`.
2. **merge** — `pmerge_sort.py` merges every chunk into
   `<out_dir>/final_sorted_data.csv` (plus a binary `.npy`).
3. **thin** — `thin_sorted_pvalues.py --thinning-factor F --method python`
   writes `<out_dir>/thinned_final_sorted_data.csv`.

Input assumptions (BGenie defaults):

- space-delimited
- header contains `chr` and `pos` columns
- one or more columns matching `col_pattern` (defaults to every `-log10p`
  column, i.e. values already on the −log10 scale)
- plain text or `.gz` — handled transparently

Restricting the set of columns:

| `col_pattern`                | thinned columns                |
| ---------------------------- | ------------------------------ |
| `-log10p$` *(default)*       | every `-log10p` column         |
| `_true-log10p$`              | only measured-trait columns    |
| `LV_.*-log10p$`              | only LV columns                |
| `(LV_.*\|.*_true)-log10p$`   | LVs + measured                 |

### Independent thinning of phenotype subgroups

By default the selected columns are partitioned into **independently-thinned
groups**, each with its own rank-1 baseline and logarithmic ladder, written to a
separate `thinned_final_sorted_data_<group>.csv`. The `groups` argument is a
comma-separated list of `name=regex` pairs; each column joins the **first** group
whose regex matches it (matched against the full column name), and a column
matching **no** group aborts the run (catching typos).

| `groups` value                          | result                                   |
| --------------------------------------- | ---------------------------------------- |
| `dTIF=_pred,mTIF=_true,LV=^LV_` *(default)* | three files: dTIFs, mTIFs, LVs        |
| `''` (empty)                            | one pooled `thinned_final_sorted_data.csv` |

The same `groups` spec must reach both stages; both `run_pipeline.sh` and the
SLURM driver pass it through automatically. `unthinned_rank` in each output is
the rank **within that group**.

Memory of the sort step scales roughly as `chunksize × n_phenotypes × 4 bytes`.
For ~1000 phenotypes, `chunksize=500000` keeps it around 2 GB resident; halve
it if you OOM.

## Running on Urblauna (left-eye May-2026 revision)

The BGenie outputs for the May-2026 left-eye revision live at:

```
/scratch/<user>/retina/GWAS/output/RunGWAS/2026_05_26_left_eye_mTIFs_dTIFs_LVs_for_revisions
```

Each chromosome file is space-delimited and carries **1058 `*-log10p` phenotype
columns**: 17 measured TIFs (`*_true-log10p`), 17 deep-TIFs (`*_pred-log10p`),
and 1024 LVs (`LV_*-log10p`). These map exactly onto the three default groups
(`mTIF` / `dTIF` / `LV`), so the default run thins each independently. To process
all of them:

```bash
# clone once on Urblauna
git clone git@github.com:BergmannLab/logarithmic-thinning.git ~/logarithmic-thinning
cd ~/logarithmic-thinning

# run the pipeline
INPUT=/scratch/<user>/retina/GWAS/output/RunGWAS/2026_05_26_left_eye_mTIFs_dTIFs_LVs_for_revisions
OUTPUT=/scratch/<user>/retina/GWAS/output/logthin_left_eye_2026_05_26

./run_pipeline.sh "$INPUT" "$OUTPUT" 1.0003 'output_ukb_imp_chr*_v3.txt' 500000
```

The default `col_pattern` picks up every `-log10p` column, so this thins all
1058 phenotypes, split into the three default groups (dTIFs, mTIFs, LVs) →
`thinned_final_sorted_data_{dTIF,mTIF,LV}.csv`. To restrict to the 17 measured
TIFs only, append a 6th argument: `'_true-log10p$'`. To thin everything together
as a single pooled output, append `''` as the 7th argument to disable grouping.

`run_pipeline.sh` runs the sort step serially across the 22 chromosomes on a
single node — fine for a quick check, but on Urblauna prefer the SLURM driver
below. Outputs (identical for both drivers):

- `$OUTPUT/chunks/` — sorted `.npz` chunks (intermediate; can be deleted after)
- `$OUTPUT/thinned_final_sorted_data_<group>.csv` — one thinned table per group
  (or a single `thinned_final_sorted_data.csv` when grouping is disabled). The
  full untrimmed ordering is never materialised on disk.

### SLURM (recommended on Urblauna)

```bash
./run_pipeline_slurm.sh "$INPUT" "$OUTPUT" 1.0003 'output_ukb_imp_chr*_v3.txt' 500000
```

Same CLI as `run_pipeline.sh`. Submits a **job array** (one task per matching
input file) for the sort step, then a merge job (`afterok` on the array), then
a thin job (`afterok` on the merge). If any sort task fails the merge and thin
are auto-cancelled — no silent partial runs.

Generated artefacts under `$OUTPUT/slurm/`:

- `inputs.txt`, `pp_columns.txt` — the per-task inputs and shared phenotype list
- `sort.sbatch`, `merge.sbatch`, `thin.sbatch` — the generated job scripts (kept
  for re-running a single stage or auditing)
- `logs/` — per-job stdout/stderr (`sort-<jobid>_<task>.{out,err}`, etc.)

Defaults: `--account=<account>`, `--partition=urblauna`,
sort 8 CPUs / 200G / 2h, merge 2 CPUs / 400G / 6h, thin 1 CPU / 100G / 1h. Edit
the `SBATCH_*` variables at the top of `run_pipeline_slurm.sh` to override.

The driver checks that all input files share the same `col_pattern` columns
before submitting, and refuses to run if `$OUTPUT/chunks/` already contains
chunks from a previous attempt.

## Output schema

`thinned_final_sorted_data.csv` columns:

- `chr` — chromosome
- `pos` — genomic position
- `pp` — `-log10(p-value)` (higher = more significant; `pp=8` ⇔ p=1e-8)
- `original_index` — row index in the merged sorted file
- `unthinned_rank` — 1-based position in the *unthinned* sorted list

`thin_sorted_pvalues_with_rows.py` is a drop-in replacement that also writes
the retained row indices to a separate file.

## Manual (per-stage) invocation

Each stage can be run on its own for tighter control, or re-run individually
by invoking the generated sbatch scripts under `$OUTPUT/slurm/`:

```bash
# 1. sort one chromosome (repeat for chr02…chr22 with the SAME --shared_dir)
python sort_merge/sort.py /path/to/chr01.txt \
  --shared_dir /tmp/sorted_chunks \
  --pp-columns trait1_true-log10p,trait2_true-log10p \
  --pp-threshold 0 \
  --chunksize 500000

# 2. merge all chunks into one sorted file
python sort_merge/pmerge_sort.py \
  --input_dir /tmp/sorted_chunks --output_dir /tmp/sort_out

# 3. apply logarithmic thinning
python thin_sorted_pvalues.py /tmp/sort_out/final_sorted_data.csv \
  --thinning-factor 1.0003 --method python \
  --output /tmp/sort_out/thinned_final_sorted_data.csv
```

Use `--p-columns` (instead of `--pp-columns`) if your inputs hold raw p-values
rather than `-log10p`.

## C++ accelerated tools

Optional faster drop-ins for very large inputs. Build once:

```bash
make -C cpp
```

Sort one file:

```bash
cpp/build/logsort --input /path/to/bgenie_output.txt \
  --output final_sorted_data.csv --threshold 0 --chunk-rows 2000000
```

Apply thinning:

```bash
cpp/build/logthin --input final_sorted_data.csv \
  --output thinned_final_sorted_data.csv --factor 1.0003
```

The C++ sorter also emits a binary stream (`final_sorted_data.bin`) for
downstream tooling. It processes a single input file, so for multi-chromosome
runs either use the Python `run_pipeline.sh` (which handles cross-file merging
through the chunks dir) or run `logsort` per chromosome and merge the outputs
yourself.
