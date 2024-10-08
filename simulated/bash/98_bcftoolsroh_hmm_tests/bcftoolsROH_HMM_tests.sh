#!/bin/bash

#SBATCH --job-name=bcf_hmm_tests
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH --mem=32000
#SBATCH --partition=general
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=avrilharder@gmail.com


##  Set username
USER=avrilh

## Set project name
PROJ=bcftoolsroh_hmm_tests

## Create a directory on /scratch
mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


# -----------------------------------------------------------------------------
# Load modules, copy data, set up variables
# -----------------------------------------------------------------------------
module load bcftools/1.17
module load samtools/1.17

mkdir data
mkdir af_files
mkdir rg_files
mkdir vt_files
mkdir output
# cp /home/amh0254/roh_param_project/files_to_archive/sample*.vcf.gz ./data/

declare -a cvgX=(05x 10x 15x 30x)
declare -a hw2az=(2.00e-6 6.7e-5 6.7e-6 6.7e-7 6.7e-8 6.7e-9 6.7e-10 6.7e-11)
declare -a az2hw=(1.16e-5 5e-4 5e-5 5e-6 5e-7 5e-8 5e-9 5e-10 5e-11 5e-12)
declare -a vit=(1e-10 1e-12)

# -----------------------------------------------------------------------------
# Generate allele frequency files for use with bcftools roh & index with tabix
# -----------------------------------------------------------------------------


# for c in ${cvgX[@]}; do
# 	
# 	bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
# 	./data/sample_pop_100_cvg_${c}_filtered_SNPs.vcf.gz \
# 	| bgzip -c > ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz
# 
# 	tabix -s1 -b2 -e2 ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz
# 
# done

# -----------------------------------------------------------------------------
# ROH calling - Genotypes and Genotype likelihoods
# +
# Extract information on ROHS i.e., exlcude information on individual
# sites
# -----------------------------------------------------------------------------
for c in ${cvgX[@]}; do
	
	# -------------------------------------------------------------------------
	# Straight defaults
	# -------------------------------------------------------------------------
	## Genotypes
	bcftools roh \
		--GTs-only 30 \
		--threads 20 \
		-o ./output/GT_sample_pop_100_cvg_${c}_defaults.txt \
		--AF-file ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz \
		./data/sample_pop_100_cvg_${c}_filtered_SNPs.vcf.gz
		
	grep "^RG" ./output/GT_sample_pop_100_cvg_${c}_defaults.txt >\
	./rg_files/GT_sample_pop_100_cvg_${c}_defaults_RG_ONLY.txt

	## Genotype likelihoods
	bcftools roh \
		--threads 20 \
		-o ./output/PL_sample_pop_100_cvg_${c}_defaults.txt \
		--AF-file ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz \
		./data/sample_pop_100_cvg_${c}_filtered_SNPs.vcf.gz
		
	grep "^RG" ./output/PL_sample_pop_100_cvg_${c}_defaults.txt >\
	./rg_files/PL_sample_pop_100_cvg_${c}_defaults_RG_ONLY.txt

	# -------------------------------------------------------------------------
	# Testing Viterbi training + options
	# -------------------------------------------------------------------------

# 	for v in ${vit[@]}; do
# 	
# 		## Genotypes
# 		bcftools roh \
# 			--GTs-only 30 \
# 			--threads 20 \
# 			-V ${v} \
# 			-o ./output/GT_sample_pop_100_cvg_${c}_viterbi_${v}.txt \
# 			--AF-file ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz \
# 			./data/sample_pop_100_cvg_${c}_filtered_SNPs.vcf.gz
# 	
# 		grep "^VT" ./output/GT_sample_pop_100_cvg_${c}_viterbi_${v}.txt >\
# 		./vt_files/GT_sample_pop_100_cvg_${c}_${v}_VT_ONLY.txt
# 			
# 	
# 		## Genotype likelihoods
# 		bcftools roh \
# 			--threads 20 \
# 			-V ${v} \
# 			-o ./output/PL_sample_pop_100_cvg_${c}_viterbi_${v}.txt \
# 			--AF-file ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz \
# 			./data/sample_pop_100_cvg_${c}_filtered_SNPs.vcf.gz
# 	
# 		grep "^VT" ./output/PL_sample_pop_100_cvg_${c}_viterbi_${v}.txt >\
# 		./vt_files/PL_sample_pop_100_cvg_${c}_${v}_VT_ONLY.txt
# 	
# 	done

	# -------------------------------------------------------------------------
	# Testing combinations of specified HMM transition probabilities
	# -------------------------------------------------------------------------

# 	for a in ${hw2az[@]}; do
# 		for h in ${az2hw[@]}; do
# 			
# 			## Genotypes
# 			bcftools roh \
# 				--GTs-only 30 \
# 				--threads 20 \
# 				-a ${a} \
# 				-H ${h} \
# 				-o ./output/GT_sample_pop_100_cvg_${c}_hw2az_${a}_az2hw_${h}.txt \
# 				--AF-file ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz \
# 				./data/sample_pop_100_cvg_${c}_filtered_SNPs.vcf.gz
# 				
# 			grep "^RG" ./output/GT_sample_pop_100_cvg_${c}_hw2az_${a}_az2hw_${h}.txt >\
# 			./rg_files/GT_sample_pop_100_cvg_${c}_hw2az_${a}_az2hw_${h}_RG_ONLY.txt
# 
# 			## Genotype likelihoods
# 			bcftools roh \
# 				--threads 20 \
# 				-a ${a} \
# 				-H ${h} \
# 				-o ./output/PL_sample_pop_100_cvg_${c}_hw2az_${a}_az2hw_${h}.txt \
# 				--AF-file ./af_files/sample_pop_100_cvg_${c}_filtered_SNPs.tab.gz \
# 				./data/sample_pop_100_cvg_${c}_filtered_SNPs.vcf.gz
# 				
# 			grep "^RG" ./output/PL_sample_pop_100_cvg_${c}_hw2az_${a}_az2hw_${h}.txt >\
# 			./rg_files/PL_sample_pop_100_cvg_${c}_hw2az_${a}_az2hw_${h}_RG_ONLY.txt
# 				
# 		done
# 	done	
done