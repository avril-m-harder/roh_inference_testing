#!/bin/bash

#SBATCH --job-name=small_05c_combine_and_genotype_gvcfs
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 07:00:00
#SBATCH --mem=16000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05d_combine_and_genotype_gvcfs

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load gatk/4.1.9.0

# -----------------------------------------------------------------------------
# Genotype GVCFs across all samples simultaneously
# -----------------------------------------------------------------------------

for c in ${cvgX[@]}; do 

	cd ${OUTPUT_DIR}

	# Set and create map file for this sample set
	MAP_FILE=${OUTPUT_DIR}/small_sample_cvg_${c}_map.list
	ls small_sample_cvg_${c}/*.g.vcf > ${MAP_FILE}
	

	OUT_FILE=small_sample_cvg_${c}_combined_gvcfs.vcf
	COMBINED_GVCFS_FILE=${OUT_FILE}

	start_logging "gatk CombineGVCFs - ${OUT_FILE}"

	gatk CombineGVCFs \
		-R ${REF_GENOME_FILE} \
		-V ${MAP_FILE} \
		-O ${OUTPUT_DIR}/${OUT_FILE}

	stop_logging

	# -------------------------------------------------------------------------
	# gatk GenotypeGVCFs
	# -------------------------------------------------------------------------

	OUT_FILE=small_sample_cvg_${c}_genotyped_gvcfs.vcf
	GENOTYPED_FILE=${OUT_FILE}

	start_logging "gatk GenotypeGVCFs- ${OUT_FILE}"

	gatk GenotypeGVCFs \
		-R ${REF_GENOME_FILE} \
		-V ${OUTPUT_DIR}/${COMBINED_GVCFS_FILE} \
		-O ${OUTPUT_DIR}/${OUT_FILE}

	stop_logging

done