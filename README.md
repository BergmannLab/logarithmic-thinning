# logarithmic-thinning

Plot a genome-wide p-value distribution without plotting every p-value.

A Manhattan or QQ plot of a large GWAS is dominated by millions of uninteresting points
that all land on top of each other near the null. **Logarithmic thinning** keeps every one
of the most significant results and then samples the rest geometrically, so the resulting
figure is visually indistinguishable from one drawn with the full data at a fraction of the
size. Thinning 171 151 566 rows at the default factor leaves 38 083.

The rule is a pure function of two numbers -- the row count and the thinning factor. No
p-values, positions or genotypes enter into it, so the selection is deterministic and
reproducible on any machine.

## Install

```bash
pip install .                 # the thinning rule alone -- no dependencies
pip install '.[pipeline]'     # + numpy/pandas, for the sort/merge pipeline
pip install '.[dev]'          # + pytest
```

The core is deliberately **dependency-free**: `logarithmic_thinning` imports only the
standard library, so anything that just needs the ranks installs nothing else.

## Use it as a library

```python
>>> from logarithmic_thinning import logarithmic_thinning, determine_threshold
>>> ranks = logarithmic_thinning(7_834_624, 1.0005)   # 1-based ranks to keep
>>> len(ranks)
17705
>>> determine_threshold(1.0005)                        # kept unthinned above this rank
2001
```

Given `total_rows` p-values **already sorted most-significant-first**, every rank above
`factor / (factor - 1)` is kept intact, and below that the ranks are sampled by repeatedly
dividing by `factor`. Index your own sorted vector with the result -- that is the whole
method.

A command-line equivalent for an existing sorted table:

```bash
thin-sorted-pvalues final_sorted_data.csv --thinning-factor 1.0003 --method python \
  --output thinned_final_sorted_data.csv
```

## Use it as a pipeline

For inputs too large to sort in memory, the repository also carries the pipeline that
produced the tables behind the published figure: sort per-chromosome GWAS output into
chunks, merge the chunks into one globally-sorted stream, and thin that stream inline.

```bash
./run_pipeline.sh <gwas_dir> <out_dir> \
                  [thinning_factor=1.0003] \
                  [glob='*'] \
                  [chunksize=500000] \
                  [col_pattern='-log10p$'] \
                  [groups='dTIF=_pred,mTIF=_true,LV=^LV_']
```

1. **sort** -- for every file matching `<gwas_dir>/<glob>`, derives the phenotype columns
   from the header (every column name matching `col_pattern`) and passes them to
   `sort_merge/sort.py`. Sorted `.npz` chunks land in `<out_dir>/chunks/`.
2. **merge + thin** -- `sort_merge/pmerge_sort.py` streams a memory-bounded k-way merge and
   applies the thinning inline, writing only the kept rows. **The full sorted ordering is
   never materialised on disk** -- that is what keeps a ~770 GB intermediate off the
   filesystem.

Input assumptions (BGenie defaults, but nothing is BGenie-specific):

- space-delimited, plain text or `.gz` (handled transparently)
- a header containing `chr` and `pos` columns
- one or more columns matching `col_pattern` -- by default every `-log10p` column, i.e.
  values already on the -log10 scale. Use `--p-columns` instead of `--pp-columns` if your
  inputs hold raw p-values.

Restricting the set of columns:

| `col_pattern`                | thinned columns                |
| ---------------------------- | ------------------------------ |
| `-log10p$` *(default)*       | every `-log10p` column         |
| `_true-log10p$`              | only measured-trait columns    |
| `LV_.*-log10p$`              | only latent-variable columns   |
| `(LV_.*\|.*_true)-log10p$`   | both of the above              |

### Independently thinned subgroups

Phenotypes measured on different scales should not share one ladder. By default the
selected columns are partitioned into **independently-thinned groups**, each with its own
rank-1 baseline, written to a separate `thinned_final_sorted_data_<group>.csv`. The
`groups` argument is a comma-separated list of `name=regex` pairs; each column joins the
**first** group whose regex matches it, and a column matching **no** group aborts the run
(which is what catches typos).

| `groups` value                              | result                                     |
| ------------------------------------------- | ------------------------------------------ |
| `dTIF=_pred,mTIF=_true,LV=^LV_` *(default)* | three files, one per group                 |
| `''` (empty)                                | one pooled `thinned_final_sorted_data.csv` |

The same spec must reach both stages; both drivers pass it through automatically.
`unthinned_rank` in each output is the rank **within that group**.

Sort-step memory scales roughly as `chunksize x n_phenotypes x 4 bytes`. For ~1000
phenotypes, `chunksize=500000` keeps it near 2 GB resident; halve it if you OOM.

### On a SLURM cluster

```bash
SBATCH_ACCOUNT=my_account SBATCH_PARTITION=my_partition \
  ./run_pipeline_slurm.sh <gwas_dir> <out_dir> [thinning_factor] [glob] [chunksize]
```

Same CLI as `run_pipeline.sh`. Submits a **job array** (one task per input file) for the
sort, then a merge job with `afterok` on the array -- if any sort task fails the merge
never starts, so there are no silent partial runs. Thinning is fused into the merge; there
is no separate thin job.

Three environment variables configure the site, all empty by default so the driver runs on
any cluster; an empty value omits the directive rather than passing an empty one:

| variable | effect |
| --- | --- |
| `SBATCH_ACCOUNT` | emits `#SBATCH --account=...` |
| `SBATCH_PARTITION` | emits `#SBATCH --partition=...` |
| `SBATCH_MODULES` | shell prologue run before `python` (module loads, venv activation, ...) |

Resource defaults -- sort 8 CPUs / 200 G / 2 h, merge 1 CPU / 96 G / 4 h -- are the
`#SBATCH` lines in the generated scripts; edit `run_pipeline_slurm.sh` to change them. The
merge's peak RAM is bounded by `--mem-budget-gb` and is independent of dataset size.

Generated under `<out_dir>/slurm/`: `inputs.txt` and `pp_columns.txt` (the per-task inputs
and shared phenotype list), `sort.sbatch` and `merge.sbatch` (kept, so a single stage can
be re-run or audited), and `logs/`. The driver checks that all inputs share the same
columns before submitting, and refuses to run if `<out_dir>/chunks/` already holds chunks
from a previous attempt.

## Output schema

`thinned_final_sorted_data*.csv`:

| column | meaning |
| --- | --- |
| `chr` | chromosome |
| `pos` | genomic position |
| `pp` | `-log10(p)` -- higher is more significant, so `pp=8` is p=1e-8 |
| `original_index` | row index in the merged sorted stream |
| `unthinned_rank` | 1-based position in the *unthinned* sorted list (within the group) |

`thin_sorted_pvalues_with_rows.py` is a drop-in replacement that also writes the retained
row indices to a separate file.

## Per-stage invocation

Each stage runs on its own for tighter control, or re-run individually through the
generated sbatch scripts:

```bash
# 1. sort one chromosome (repeat, with the SAME --shared_dir)
python sort_merge/sort.py /path/to/chr01.txt \
  --shared_dir /tmp/sorted_chunks \
  --pp-columns trait1_true-log10p,trait2_true-log10p \
  --pp-threshold 0 \
  --chunksize 500000

# 2. merge all chunks, thinning inline
python sort_merge/pmerge_sort.py \
  --input_dir /tmp/sorted_chunks --output_dir /tmp/sort_out \
  --thinning-factor 1.0003

# 3. or thin an already-sorted table on its own
python thin_sorted_pvalues.py /tmp/sort_out/final_sorted_data.csv \
  --thinning-factor 1.0003 --method python \
  --output /tmp/sort_out/thinned_final_sorted_data.csv
```

## C++ accelerated tools

Optional faster drop-ins for very large inputs. Build once with `make -C cpp`:

```bash
cpp/build/logsort --input /path/to/gwas_output.txt \
  --output final_sorted_data.csv --threshold 0 --chunk-rows 2000000

cpp/build/logthin --input final_sorted_data.csv \
  --output thinned_final_sorted_data.csv --factor 1.0003
```

The C++ sorter also emits a binary stream (`final_sorted_data.bin`). It processes a single
input file, so for multi-chromosome runs either use `run_pipeline.sh` (which handles
cross-file merging through the chunks dir) or run `logsort` per chromosome and merge
yourself.

## Tests

```bash
pip install '.[dev]' && pytest
```

- `test_published_figure_inputs.py` -- pins the rank selection behind the published figure
  by checksum. Runs anywhere, needs no data.
- `test_minimal_pipeline.py`, `test_group_thinning.py` -- the pipeline and the grouping on
  small fixtures.
- `test_large_pipeline_consistency.py` -- the Python and C++ pipelines agree on a synthetic
  genome with planted peaks.
- `test_cpp_tools.py` -- the C++ tools alone (builds `cpp/` on demand).

Every test file is also runnable as a plain script.

## Provenance

This is the thinning method used for the Manhattan panel of a retinal-vasculometry GWAS,
where the three trait sets were thinned at factor 1.0003 from row counts 2 832 855,
2 914 788 and 171 151 566. `examples/reproduce_figure3a_ranks.py` checks a set of
published thinned tables against a freshly computed selection.

Note the scope of that check, and of the test above: both verify **which rows** were kept,
not the contents of those rows. Regenerating the contents means running the pipeline end to
end. The two are different claims and only the first is made here.

## Licence

MIT -- see [LICENSE](LICENSE). Deliberately permissive: this is a small general-purpose
utility, and MIT keeps it usable regardless of the licence of whatever it is plotted into.
