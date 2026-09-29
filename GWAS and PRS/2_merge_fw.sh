#!/bin/bash
#SBATCH --job-name=gwas             # Job name
#SBATCH --output=%j.out             # Stdout (%j expands to jobId)
#SBATCH --error=%j.err              # Stderr (%j expands to jobId)
#SBATCH --partition=share
#SBATCH --nodes=1                   # Number of nodes requested
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00             # walltime

ukbTrait=FW
plinktest=glm.*

cat ukb_chr1.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'BEGIN {print "SNP", "CHR", "BP", "Ref", "Alt", "A1", "EAF", "N", "Beta", "Se", "P", "Err"}; NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >  ukb_${ukbTrait}.assoc
for chr in {2..22}; do
cat ukb_chr${chr}.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >>  ukb_${ukbTrait}.assoc
done


ukbTrait=FW_Projection
plinktest=glm.*

cat ukb_chr1.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'BEGIN {print "SNP", "CHR", "BP", "Ref", "Alt", "A1", "EAF", "N", "Beta", "Se", "P", "Err"}; NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >  ukb_${ukbTrait}.assoc
for chr in {2..22}; do
cat ukb_chr${chr}.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >>  ukb_${ukbTrait}.assoc
done

ukbTrait=FW_Commissural
plinktest=glm.*

cat ukb_chr1.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'BEGIN {print "SNP", "CHR", "BP", "Ref", "Alt", "A1", "EAF", "N", "Beta", "Se", "P", "Err"}; NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >  ukb_${ukbTrait}.assoc
for chr in {2..22}; do
cat ukb_chr${chr}.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >>  ukb_${ukbTrait}.assoc
done

ukbTrait=FW_Association
plinktest=glm.*

cat ukb_chr1.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'BEGIN {print "SNP", "CHR", "BP", "Ref", "Alt", "A1", "EAF", "N", "Beta", "Se", "P", "Err"}; NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >  ukb_${ukbTrait}.assoc
for chr in {2..22}; do
cat ukb_chr${chr}.${ukbTrait}.${plinktest} | awk -v FS="\t" -v OFS="\t" 'NR > 1  {print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >>  ukb_${ukbTrait}.assoc
done
