# ADMIXTURE 

## Introduction
This page is a work in progress!
This repo explains how to run ADMIXTURE on Linux and R from a .vcf file. This pipeline runs in Linux and R and relies on BCFtools (version BCFtools/1.21-GCC-13.3.0), PLINK (version PLINK/2.00a3.7-gfbf-2023a), and VCFtools (version VCFtools/0.1.16-GCC-13.2.0). For this entire pipeline, it may be necessary to split each step up into smaller sample groups. All code should be run in scratch.

Contact: Camille Block (camilleblock@vt.edu)

## Remember to create a SLURM Script
```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --job-name STACKS
#SBATCH --account=bedbug
#SBATCH --partition=normal_q
```
## Load needed directories 
```bash
module load BCFtools/1.21-GCC-13.3.0
module load PLINK/2.00a3.7-gfbf-2023a
module load VCFtools/0.1.16-GCC-13.2.0
```
## Subset vcf file so that it runs faster (break up samples into groups of 30-50)
```bash
bcftools view -S keep1.txt snps_merged.integer.maf001.geno75.final.vcf.gz -Oz -o subset1.vcf.gz
```
## Index files so that Plink can read them (this must be done separately for each group of samples)
```bash
tabix -p vcf subset1.vcf.gz
```
## Convert each vcf to bed files (tab delimited file that stores chromosomes as coordinates) and filter 
```bash
plink --bfile snps1 \
     --maf 0.001 \
     --geno 0.75 \
     --allow-extra-chr \
     --make-bed \
     --out snps_subset1.int75
```
## Load ADMIXTURE 
```bash
export PATH=/home/camilleblock/.local/easybuild/software/ADMIXTURE/1.3.0-x86_64:$PATH
echo 'export PATH=/home/camilleblock/.local/easybuild/software/ADMIXTURE/1.3.0-x86_64:$PATH' >> ~/.bashrc
source ~/.bashrc
```
## Run Admixture for each bed file and for K=1 through K=8 (or more if needed)
```bash
admixture --cv snps_subset1.int75.bed 1 > logsubset1_1.out
```
