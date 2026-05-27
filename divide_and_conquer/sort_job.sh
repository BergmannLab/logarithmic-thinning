#!/bin/bash
#SBATCH --job-name=psort
#SBATCH --output=slurm_umapsort_chr$1-%j.out
#SBATCH --error=slurm_umapsort_chr$1-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=100G
#SBATCH --time=00:20:00
source /dcsrsoft/spack/bin/setup_dcsrsoft
module load gcc python/3.8

#INPUT_FOLDER="/scratch/<user>/retina/GWAS/output/RunGWAS/2024_09_04_RETFound_LV_Extraction_Oxford_Method_all_images_raw_weights_rebuttal_SE/"
INPUT_FOLDER="/scratch/<user>/retina/GWAS/output/RunGWAS/2024_10_2_all_pred_vs_true_rebuttal_SE/"
#INPUT_FOLDER="/scratch/<user>/retina/GWAS/output/RunGWAS/2024_11_06_lv_pcs17"
INPUT_FILE="$INPUT_FOLDER/output_ukb_imp_chr$1_v3.txt"
SHARED_DIR="$INPUT_FOLDER/true_sorted_chunks_pp2"

mkdir -p "$SHARED_DIR"

python sort.py $INPUT_FILE --shared_dir $SHARED_DIR --chunksize 5000000000 --pp_threshold 2

