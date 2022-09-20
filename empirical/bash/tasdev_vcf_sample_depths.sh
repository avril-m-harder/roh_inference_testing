#!/bin/bash

#SBATCH --job-name=depth
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

vcftools/0.1.17

cd /scratch/avrilh/rohparam_03_gatk

vcftools --vcf sorted_tasdev_all_chroms_15samps_fullcovg.recode.vcf --depth --out full
vcftools --vcf sorted_tasdev_all_chroms_15samps_5X.recode.vcf --depth --out 5X
vcftools --vcf sorted_tasdev_all_chroms_15samps_10X.recode.vcf --depth --out 10X
vcftools --vcf sorted_tasdev_all_chroms_15samps_15X.recode.vcf --depth --out 15X
vcftools --vcf sorted_tasdev_all_chroms_15samps_30X.recode.vcf --depth --out 30X
