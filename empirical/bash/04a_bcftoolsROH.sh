#!/bin/bash

#SBATCH --job-name=roh_03b_gatk
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=16000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std


##  Set username
USER=avrilh

## Set project name
PROJ=rohparam_04a_bcftoolsROH

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into directory
cd /scratch/${USER}/${PROJ}/


## --------------------------------
## Load modules
module load bcftools/1.11
module load samtools/1.11
module load vcftools/0.1.17
module load htslib/1.11


## --------------------------------
## Loop over full and downsampled VCF files
while read -a line
do

# bcftools sort \
# -o ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${line[0]}.recode.vcf \
# ../rohparam_03_gatk/final_filt_tasdev_all_chroms_15samps_${line[0]}.recode.vcf
# 
# 
# bgzip -@ 20 -c \
# ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${line[0]}.recode.vcf \
# > ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${line[0]}.recode.vcf.gz
# 
# tabix -p vcf --csi \
# ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${line[0]}.recode.vcf.gz

## --------------------------------
## Generate allele frequency files for use with bcftools roh & index
# bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
# ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${line[0]}.recode.vcf.gz \
# | bgzip -c > final_tasdev_all_chroms_15samps_${line[0]}.tab.gz
# 
# tabix -s1 -b2 -e2 --csi final_tasdev_all_chroms_15samps_${line[0]}.tab.gz

## --------------------------------
## Call ROHs with GT and PL 
# bcftools roh \
# --GTs-only 30 \
# --threads 20 \
# -o tasdev_GTonly_${line[0]}.txt \
# --AF-file final_tasdev_all_chroms_15samps_${line[0]}.tab.gz \
# ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${line[0]}.recode.vcf.gz
# 
# bcftools roh \
# --threads 20 \
# -o tasdev_PLonly_${line[0]}.txt \
# --AF-file final_tasdev_all_chroms_15samps_${line[0]}.tab.gz \
# ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${line[0]}.recode.vcf.gz


## --------------------------------
## Extract information on ROHs (i.e., exclude information on individual sites)
grep "^RG" tasdev_GTonly_${line[0]}.txt > \
tasdev_GTonly_${line[0]}_RG_ONLY.txt

grep "^RG" tasdev_PLonly_${line[0]}.txt > \
tasdev_PLonly_${line[0]}L_RG_ONLY.txt


done < /home/amh0254/roh_param_project/sample_lists/covg_levels.txt

