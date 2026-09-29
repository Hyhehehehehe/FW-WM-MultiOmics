#!/bin/bash
#SBATCH --job-name=gwas           
#SBATCH --output=%j.out           
#SBATCH --error=%j.err            
#SBATCH --partition=share
#SBATCH --nodes=1                   
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00            

chr=$1
plink2 \
--pfile ukb22828_eur_chr${chr} \
--memory 2000 \
--remove ukb_eur.rmid \
--glm hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,beta,se,p,err \
--pheno FW.txt \
--covar FW.covar \
--covar-variance-standardize \
--covar-name Age_attending_2 sex Mean_SBP_2 gen_batch gen_PC1 gen_PC2 gen_PC3 gen_PC4 gen_PC5 gen_PC6 gen_PC7 gen_PC8 gen_PC9 gen_PC10 \
--extract ukb_eur5k.bim \
--out ukb_chr${chr}


