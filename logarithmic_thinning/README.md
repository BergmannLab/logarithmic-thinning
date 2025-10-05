
# Logarithmic Thinning for the RETFound Project

This pipeline applies logarithmic thinning to GWAS p-values, starting from BGenie output files. The process is optimized for SLURM clusters and is designed to efficiently reduce the number of p-values for downstream analysis, while keeping the most significant results and a representative sample of less significant ones.

## Steps

1. **Sort p-values, keeping only those below a certain threshold (to avoid excessive data size and speed up processing)**
	 - In `divide_and_conquer/sort_job.sh`, set the BGenie output folder (`INPUT_FOLDER`).
	 - Optionally, for mTIFs and dTIFs, you can subset columns to pick only those ending with `pred-log10p`.
	 - Run the sorting for all chromosomes:
		 ```bash
		 for CHROM in {1..22}; do sbatch sort_job.sh $CHROM & done
		 ```

2. **Merge-sort the chromosome-wise lists**
	 - After sorting, merge all chromosome results into a single sorted file:
		 ```bash
		 python divide_and_conquer/pmerge_sort.py --input_dir /path/to/true_sorted_chunks_pp2 --output_dir .
		 ```

3. **Apply logarithmic thinning**
	 - From the root project directory, run:
		 ```bash
		 python thin_sorted_pvalues.py divide_and_conquer/final_sorted_data.csv 1.0003 python
		 ```
	 - Output will appear in the `divide_and_conquer` sub-directory. Rename the output file for uniqueness if needed.

## Final Output

The final output is a CSV file (e.g., `thinned_python_final_sorted_data.csv`) containing the thinned, sorted p-values. This file is suitable for downstream analysis and visualization (e.g., Manhattan plots).

### Output Columns

- `chr`: Chromosome number
- `pos`: Genomic position
- `pp`: The -log10(p-value) (higher values are more significant)
- `original_index`: The original row index from the merged file
- `row_number`: The row number in the final thinned file

#### What is `pp`?
`pp` stands for the -log10(p-value) as output by BGenie. A higher `pp` means a more significant association (e.g., `pp=8` means p-value = 1e-8).

---
This workflow ensures you keep the most significant GWAS hits and a representative sample of less significant ones, while making the data manageable for further analysis.
