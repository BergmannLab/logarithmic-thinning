
# Logarithmic Thinning for the RETFound Project

This pipeline applies logarithmic thinning to GWAS p-values, starting from BGenie output files. The process is optimized for SLURM clusters and is designed to efficiently reduce the number of p-values for downstream analysis, while keeping the most significant results and a representative sample of less significant ones.

## Steps

1. **Sort the per-chromosome files, keeping only rows above a p-value threshold**
    - Configure `divide_and_conquer/sort_job.sh` with the input folder (`INPUT_FOLDER`), desired chunk size, delimiter, and the columns you want to keep.  
      You can now identify the -log10 columns explicitly via `--pp-columns` or let the sorter convert raw p-values with `--p-columns`.
    - Example (single chromosome, interactive run):
      ```bash
      python divide_and_conquer/sort.py /path/to/chr01.txt \
        --shared_dir /tmp/chr01_chunks \
        --pp-columns trait_true-log10p \
        --pp-threshold 0 \
        --chunksize 5000000
      ```
    - To scan a batch on SLURM:
      ```bash
      for CHROM in {1..22}; do
        sbatch sort_job.sh $CHROM
      done
      ```
      where `sort_job.sh` forwards the relevant CLI arguments.

2. **Merge-sort the chromosome-wise lists**
	 - After sorting, merge all chromosome results into a single sorted file:
		 ```bash
		 python divide_and_conquer/pmerge_sort.py --input_dir /path/to/true_sorted_chunks_pp2 --output_dir .
		 ```

3. **Apply logarithmic thinning**
	 - From the root project directory, run:
		 ```bash
		 python thin_sorted_pvalues.py divide_and_conquer/final_sorted_data.csv \
		   --thinning-factor 1.0003 --method python
		 ```
	 - If you also want the retained row indices in a separate file:
		 ```bash
		 python thin_sorted_pvalues_with_rows.py divide_and_conquer/final_sorted_data.csv \
		   --thinning-factor 1.0003 --method python
		 ```
	 - Output will appear in the `divide_and_conquer` sub-directory. Rename the output file for uniqueness if needed.

## C++ Accelerated Tools

For larger workloads you can use the C++ implementations that mirror the Python pipeline:

- Build the executables once:
  ```bash
  make -C cpp
  ```
- Run the sorter (chunk rows defaults to 1,000,000):
  ```bash
  cpp/build/logsort --input path/to/bgenie_output.txt --output final_sorted_data.csv --threshold 0 --chunk-rows 2000000
  ```
  Temporary chunk files are written next to the output unless `--tmpdir` is provided.
- Apply thinning on the sorted CSV:
  ```bash
  cpp/build/logthin --input final_sorted_data.csv --output thinned_final_sorted_data.csv --factor 1.0003
  ```

The C++ sorter also emits a binary stream (`final_sorted_data.bin`) that contains the sorted records in the order written to the CSV for downstream tooling.

## Final Output

The final output is a CSV file (e.g., `thinned_python_final_sorted_data.csv`) containing the thinned, sorted p-values. This file is suitable for downstream analysis and visualization (e.g., Manhattan plots).

### Output Columns

- `chr`: Chromosome number
- `pos`: Genomic position
- `pp`: The -log10(p-value) (higher values are more significant)
- `original_index`: The original row index from the merged file
- `unthinned_rank`: Position within the *unthinned* sorted list (1-based)

#### What is `pp`?
`pp` stands for the -log10(p-value) as output by BGenie. A higher `pp` means a more significant association (e.g., `pp=8` means p-value = 1e-8).

---
This workflow ensures you keep the most significant GWAS hits and a representative sample of less significant ones, while making the data manageable for further analysis.
