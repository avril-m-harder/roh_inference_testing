#!/bin/bash

#SBATCH --job-name=_demo__05e_filter_variants
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 30:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05e_filter_variants

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load gatk/4.1.9.0


# ------------------------------------------------------------------------------
# Flag Variants that do not pass filters
# ------------------------------------------------------------------------------

for c in ${cvgX[@]}; do 

	OUT_FILE=_demo__sample_cvg_${c}_flagged_gvcfs.vcf
	FLAGGED_FILE=${OUT_FILE}
	GENOTYPED_FILE=_demo__sample_cvg_${c}_genotyped_gvcfs.vcf

	start_logging "gatk VariantFiltration- ${OUT_FILE}"

	gatk VariantFiltration \
		-R ${REF_GENOME_FILE} \
		-V ${OUTPUT_DIR}/${GENOTYPED_FILE} \
		-O ${OUTPUT_DIR}/${OUT_FILE} \
		--filter-name "QD" \
		--filter-expression "QD < 2.0" \
		--filter-name "FS" \
		--filter-expression "FS > 40.0" \
		--filter-name "SOR" \
		--filter-expression "SOR > 5.0" \
		--filter-name "MQ" \
		--filter-expression "MQ < 20.0" \
		--filter-name "MQRankSum" \
		--filter-expression " MQRankSum < -3.0 || MQRankSum > 3.0" \
		--filter-name "ReadPosRankSum" \
		--filter-expression "ReadPosRankSum < -8.0"

	stop_logging

	# ------------------------------------------------------------------------------
	# Keep only SNPs that pass the above filters
	# ------------------------------------------------------------------------------

	OUT_FILE=_demo__sample_cvg_${c}_filtered_SNPS.vcf

	start_logging "gatk SelectVariants - ${OUT_FILE}"

	gatk SelectVariants \
		-R ${REF_GENOME_FILE} \
		-V ${OUTPUT_DIR}/${FLAGGED_FILE} \
		--select-type-to-include SNP \
		-select 'vc.isNotFiltered()' \
		-O ${OUTPUT_DIR}/${OUT_FILE}

	stop_logging
	
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

cp _demo__sample_cvg_*_filtered_SNPS.vcf \
/home/amh0254/roh_param_project/roh_inference_testing/simulated/data/05_var_call/
