#!/bin/bash

#SBATCH --job-name=roh_04_bcftools
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=16000
#SBATCH -t 300:00:00
#SBATCH --output=04_bcftools-%j.out
#SBATCH --error=error-04_bcftools-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std


##  Set username
USER=avrilh

## Set project name
PROJ=archive_first_roh_param_runs/rohparam_04a_bcftoolsROH_empirical

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into directory
cd /scratch/${USER}/${PROJ}/


# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load bcftools/1.17
module load samtools/1.17
module load vcftools/0.1.17
module load htslib/1.17
module load R/4.2.3

declare -a cvgs=("5X" "10X" "15X" "30X" "fullcovg")
declare -a methods=("PL" "GT")


# -----------------------------------------------------------------------------
# Prep VCF files
# -----------------------------------------------------------------------------

# for c in ${cvgs[@]}
# do
# 
# 	bcftools sort \
# 	-o ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${c}.recode.vcf \
# 	../rohparam_03_gatk/final_filt_tasdev_all_chroms_15samps_${c}.recode.vcf
# 	
# 	
# 	bgzip -@ 20 -c \
# 	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${c}.recode.vcf \
# 	> ../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${c}.recode.vcf.gz
# 	
# 	tabix -p vcf --csi \
# 	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${c}.recode.vcf.gz


	# -----------------------------------------------------------------------------
	# Generate allele frequency files for use with bcftools roh & index
	# -----------------------------------------------------------------------------

# 	bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
# 	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${c}.recode.vcf.gz \
# 	| bgzip -c > final_tasdev_all_chroms_15samps_${c}.tab.gz
# 
# 	tabix -s1 -b2 -e2 --csi final_tasdev_all_chroms_15samps_${c}.tab.gz
# 
# 
# 	# -----------------------------------------------------------------------------
# 	# Run Viterbi training to determine appropriate HMM transition probabilities 
# 	# -----------------------------------------------------------------------------
# 
# 	## Genotypes
# 	bcftools roh \
# 	--GTs-only 30 \
# 	--threads 20 \
# 	--viterbi-training 1e-10 \
# 	-o tasdev_GTonly_${c}_viterbi.txt \
# 	--AF-file final_tasdev_all_chroms_15samps_${c}.tab.gz \
# 	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${c}.recode.vcf.gz
# 
# 	## Genotype likelihoods
# 	bcftools roh \
# 	--threads 20 \
# 	--viterbi-training 1e-10 \
# 	-o tasdev_PLonly_${c}_viterbi.txt \
# 	--AF-file final_tasdev_all_chroms_15samps_${c}.tab.gz \
# 	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${c}.recode.vcf.gz
# 
# 	grep "^VT" tasdev_GTonly_${c}_viterbi.txt > \
# 	tasdev_GTonly_${c}_VT_ONLY.txt
# 
# 	grep "^VT" tasdev_PLonly_${c}_viterbi.txt > \
# 	tasdev_PLonly_${c}_VT_ONLY.txt
# 
# done


# -----------------------------------------------------------------------------
# Use Viterbi training output to calculate HMM transition probabilities 
# -----------------------------------------------------------------------------

Rscript \
/home/amh0254/roh_param_project/roh_inference_testing/empirical/R/\
EASLEY_01_evaluate_viterbi_estimates.R


# -----------------------------------------------------------------------------
# Call ROHs with Viterbi-trained and default transition probabilities
# -----------------------------------------------------------------------------
while read -a line
do
	COVG=${line[0]}
	HW2AZ=${line[1]}
	AZ2HW=${line[3]}
	
	## defaults
	bcftools roh \
	--GTs-only 30 \
	--threads 20 \
	-o tasdev_GTonly_${COVG}_default.txt \
	--AF-file final_tasdev_all_chroms_15samps_${COVG}.tab.gz \
	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${COVG}.recode.vcf.gz
	
	## Viterbi-trained values
	bcftools roh \
	--GTs-only 30 \
	--threads 20 \
	--hw-to-az ${HW2AZ} \
	--az-to-hw ${AZ2HW} \
	-o tasdev_GTonly_${COVG}_vtrained.txt \
	--AF-file final_tasdev_all_chroms_15samps_${COVG}.tab.gz \
	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${COVG}.recode.vcf.gz
	
done < GT_viterbi_estimates.txt

while read -a line
do
	COVG=${line[0]}
	HW2AZ=${line[1]}
	AZ2HW=${line[3]}
	
	## defaults
	bcftools roh \
	--threads 20 \
	-o tasdev_PLonly_${COVG}_default.txt \
	--AF-file final_tasdev_all_chroms_15samps_${COVG}.tab.gz \
	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${COVG}.recode.vcf.gz
	
	## Viterbi-trained values
	bcftools roh \
	--threads 20 \
	--hw-to-az ${HW2AZ} \
	--az-to-hw ${AZ2HW} \
	-o tasdev_PLonly_${COVG}_vtrained.txt \
	--AF-file final_tasdev_all_chroms_15samps_${COVG}.tab.gz \
	../rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_${COVG}.recode.vcf.gz
	
done < PL_viterbi_estimates.txt


# -----------------------------------------------------------------------------
# Extract information on ROHs (i.e., exclude information on individual sites)
# -----------------------------------------------------------------------------

for m in ${methods[@]}
do
	for c in ${cvgs[@]}
	do

		grep "^RG" tasdev_${m}only_${c}_vtrained.txt > \
		tasdev_${m}only_${c}_vtrained_RG_ONLY.txt
	
		grep "^RG" tasdev_${m}only_${c}_default.txt > \
		tasdev_${m}only_${c}_default_RG_ONLY.txt
	
		
	done
done

cd /home/amh0254/roh_param_project/roh_inference_testing/empirical/bash/
mv *.out ./script_stdouts/
mv *.err ./script_stdouts/