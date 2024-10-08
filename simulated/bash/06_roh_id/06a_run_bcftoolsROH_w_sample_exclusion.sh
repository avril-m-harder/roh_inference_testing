#!/bin/bash

#SBATCH --job-name=decline_06a_run_bcftoolsROH
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH --mem=32000


# -----------------------------------------------------------------------------
# Set Variables
# -----------------------------------------------------------------------------

STEP=06_roh_id
PREV_STEP=05_var_call
SCRIPT=06a_run_bcftoolsROH

## Load variables and functions from settings file
source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

declare -a methods=("PL" "GT")


# -----------------------------------------------------------------------------
## Load modules
# -----------------------------------------------------------------------------

module load bcftools/1.17
module load samtools/1.17
module load vcftools/0.1.17
module load htslib/1.17
module load R/4.2.3


# -----------------------------------------------------------------------------
# Process final filtered SNPs
# -----------------------------------------------------------------------------

for c in ${cvgX[@]}; do

	cd ${OUTPUT_DIR} 
	
	# ----------------------------------------------------------------------
	# Generate allele frequency files for use with bcftools roh & index
	# with tabix
	# ----------------------------------------------------------------------

	IN_FILE=decline_sample_cvg_${c}_filtered_SNPs.vcf
	OUT_FILE=decline_sample_cvg_${c}_bcftools_filteredSNPs.tab.gz

	start_logging "bcftools query/tabix  - ${OUT_FILE}"
	
	## Remove 2 samples with extremely high homozygosity that hang the Viterbi step
	vcftools \
	--vcf ${INPUT_DIR}/${IN_FILE} \
	--recode \
	--recode-INFO-all \
	--remove-indv i33 \
	--remove-indv i48 \
	--out ${INPUT_DIR}/excluded_samps_decline_cvg_${c}
	
	EX_IN_FILE=excluded_samps_decline_cvg_${c}.recode.vcf
	EX_OUT_FILE=excluded_samps_decline_cvg_${c}.tab.gz

	## Extract allele frequencies from sample-excluded vcf and full vcf
	bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
	${INPUT_DIR}/${EX_IN_FILE} | bgzip -c >${OUTPUT_DIR}/${EX_OUT_FILE}

	tabix -s1 -b2 -e2 ${OUTPUT_DIR}/${EX_OUT_FILE}
	
	bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
	${INPUT_DIR}/${IN_FILE} | bgzip -c >${OUTPUT_DIR}/${OUT_FILE}

	tabix -s1 -b2 -e2 ${OUTPUT_DIR}/${OUT_FILE}

	stop_logging

		
 	# -----------------------------------------------------------------------------
 	# Run Viterbi training to determine appropriate HMM transition probabilities 
 	# -----------------------------------------------------------------------------

	## Genotypes
	IN_FILE=excluded_samps_decline_cvg_${c}.recode.vcf
	GT_OUT_FILE=decline_sample_cvg_${c}_bcftools_roh_gt_viterbi.txt
	AF_FILE=excluded_samps_decline_cvg_${c}.tab.gz

	start_logging "bcftools roh gt - ${OUT_FILE}"

	bcftools roh \
		--GTs-only 30 \
		--threads 20 \
		--viterbi-training 1e-10 \
		-o ${OUTPUT_DIR}/${GT_OUT_FILE} \
		--AF-file ${OUTPUT_DIR}/${AF_FILE} \
		${INPUT_DIR}/${IN_FILE}
	
	## Genotype likelihoods	
	IN_FILE=excluded_samps_decline_cvg_${c}.recode.vcf
	PL_OUT_FILE=decline_sample_cvg_${c}_bcftools_roh_pl_viterbi.txt
	AF_FILE=excluded_samps_decline_cvg_${c}.tab.gz

	start_logging "bcftools roh pl - ${OUT_FILE}"

	bcftools roh \
		--threads 20 \
		--viterbi-training 1e-10 \
		-o ${OUTPUT_DIR}/${PL_OUT_FILE} \
		--AF-file ${OUTPUT_DIR}/${AF_FILE} \
		${INPUT_DIR}/${IN_FILE}
		
	grep "^VT" ${OUTPUT_DIR}/${GT_OUT_FILE} > \
	${OUTPUT_DIR}/decline_sample_cvg_${c}_GT_VT_ONLY.txt

	grep "^VT" ${OUTPUT_DIR}/${PL_OUT_FILE} > \
	${OUTPUT_DIR}/decline_sample_cvg_${c}_PL_VT_ONLY.txt

done

# -----------------------------------------------------------------------------
# Use Viterbi training output to calculate HMM transition probabilities 
# -----------------------------------------------------------------------------

Rscript /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/06_roh_id/decline/EASLEY_01_evaluate_viterbi_estimates_w_make.R


# -----------------------------------------------------------------------------
# Call ROHs with Viterbi-trained and default transition probabilities
# -----------------------------------------------------------------------------
## Genotypes
while read -a line
do
	COVG=${line[0]}
	HW2AZ=${line[1]}
	AZ2HW=${line[3]}
	
	AF_FILE=decline_sample_cvg_${COVG}_bcftools_filteredSNPs.tab.gz
	IN_FILE=decline_sample_cvg_${COVG}_filtered_SNPs.vcf
	VIT_OUT_FILE=decline_sample_cvg_${COVG}_bcftools_roh_GT_vtrained.txt
	DEF_OUT_FILE=decline_sample_cvg_${COVG}_bcftools_roh_GT_defaults.txt

	## defaults
	bcftools roh \
	--GTs-only 30 \
	--threads 20 \
	-o ${OUTPUT_DIR}/${DEF_OUT_FILE} \
	--AF-file ${AF_FILE} \
	${INPUT_DIR}/${IN_FILE}

	## Viterbi-trained values
	bcftools roh \
	--GTs-only 30 \
	--threads 20 \
	--hw-to-az ${HW2AZ} \
	--az-to-hw ${AZ2HW} \
	-o ${OUTPUT_DIR}/${VIT_OUT_FILE} \
	--AF-file ${AF_FILE} \
	${INPUT_DIR}/${IN_FILE}
	
	## Extract ROH information	
	grep "^RG" ${OUTPUT_DIR}/decline_sample_cvg_${COVG}_bcftools_roh_GT_vtrained.txt > \
	${OUTPUT_DIR}/decline_GT_${COVG}_vtrained_RG_ONLY.txt
	
	grep "^RG" ${OUTPUT_DIR}/decline_sample_cvg_${COVG}_bcftools_roh_GT_defaults.txt > \
	${OUTPUT_DIR}/decline_GT_${COVG}_defaults_RG_ONLY.txt

done < decline_GT_viterbi_estimates.txt

while read -a line
do
	COVG=${line[0]}
	HW2AZ=${line[1]}
	AZ2HW=${line[3]}
	
	AF_FILE=decline_sample_cvg_${COVG}_bcftools_filteredSNPs.tab.gz
	IN_FILE=decline_sample_cvg_${COVG}_filtered_SNPs.vcf
	VIT_OUT_FILE=decline_sample_cvg_${COVG}_bcftools_roh_PL_vtrained.txt
	DEF_OUT_FILE=decline_sample_cvg_${COVG}_bcftools_roh_PL_defaults.txt

	## defaults
	bcftools roh \
	--threads 20 \
	-o ${OUTPUT_DIR}/${DEF_OUT_FILE} \
	--AF-file ${AF_FILE} \
	${INPUT_DIR}/${IN_FILE}

	## Viterbi-trained values
	bcftools roh \
	--threads 20 \
	--hw-to-az ${HW2AZ} \
	--az-to-hw ${AZ2HW} \
	-o ${OUTPUT_DIR}/${VIT_OUT_FILE} \
	--AF-file ${AF_FILE} \
	${INPUT_DIR}/${IN_FILE}
	
	## Extract ROH information	
	grep "^RG" ${OUTPUT_DIR}/decline_sample_cvg_${COVG}_bcftools_roh_PL_vtrained.txt > \
	${OUTPUT_DIR}/decline_PL_${COVG}_vtrained_RG_ONLY.txt
	
	grep "^RG" ${OUTPUT_DIR}/decline_sample_cvg_${COVG}_bcftools_roh_PL_defaults.txt > \
	${OUTPUT_DIR}/decline_PL_${COVG}_defaults_RG_ONLY.txt

done < decline_PL_viterbi_estimates.txt

	
# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

cp decline*RG_ONLY.txt \
/home/amh0254/roh_param_project/roh_inference_testing/simulated/data/06_roh_id