#!/bin/bash
#SBATCH --partition=cu,fat
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-22
#SBATCH --output="logs/%x.%A.%a.%j.out"
#SBATCH --error="logs/%x.%A.%a.%j.err"

set -euo pipefail
mkdir -p effect

export PATH="/home/xcao/miniconda3/envs/prscs/bin:$PATH"

N_THREADS=4
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

pheno="$1"
chr=${SLURM_ARRAY_TASK_ID}

ss="/data1/projects/haomicslab/xcao/cooperation/FW/PRS/fmt/ukb_${pheno}.ss"
geno="/data1/projects/haomicslab/xcao/UKB/genotype/impt/ukb22828_eur_chr${chr}"
outdir="/data1/projects/haomicslab/xcao/cooperation/FW/PRS/effect/ukb_${pheno}"

python3 /home/xcao/PRScs/PRScs.py \
        --ref_dir=/home/xcao/PRScs/ref_panel/1kg/ldblk_1kg_eur \
        --bim_prefix="${geno}" \
        --sst_file="${ss}" \
        --n_gwas=31644 \
        --out_dir="${outdir}" \
        --seed=20250529 \
        --chrom=${chr}

#for pheno in "FW_Association" "FW_Commissural" "FW_Projection" "FW"; do
#  sbatch 2.PRScs.sh ${pheno}
#done