#!/bin/bash

#SBATCH --job-name=06b_run_PLINK
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=06_roh_id
PREV_STEP=05_var_call
SCRIPT=06b_run_PLINK

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

declare -a dems=(decline large-1000)


# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load plink/1.9


# ----------------------------------------------------------------------
# Loop over coverage levels within demo scenario
# ----------------------------------------------------------------------

for d in ${dems[@]}; do 
	for c in ${cvgX[@]}; do

		cd ${OUTPUT_DIR} 

		# ----------------------------------------------------------------------
		# Convert VCF to PLINK format(s)
		# ----------------------------------------------------------------------

		IN_FILE=${d}_sample_cvg_${c}_filtered_SNPs.vcf
		OUT_FILE=${d}_sample_cvg_${c}_plink

		start_logging "Convert VCF to PLINK"

		plink \
			--vcf ${INPUT_DIR}/${IN_FILE} \
			--allow-extra-chr \
			--out ${OUTPUT_DIR}/${OUT_FILE}

		stop_logging

		# -----------------------------------------------------------------------------
		# PLINK parameters
		# -----------------------------------------------------------------------------

		declare -a phwh=(0 1 2)         # Values for -homozyg-window-het
		declare -a phwm=(2 5 50)        # Values for -homozyg-window-missing
		declare -a phws=(50 100 1000)   # Values for -homozyg-window-snp
		declare -a phzd=(50)            # Values for -homozyg-density
		declare -a phzg=(500 1000)      # Values for -homozyg-gap
		declare -a phwt=(0.01 0.05 0.1) # Values for -homozyg-window-threshold
		declare -a phzs=(10 100 1000)   # Values for -homozyg-snp
		declare -a phzk=(100)           # Values for -homozyg-kb

		# ----------------------------------------------------------------------
		# Run PLINK to ID ROHs
		# ----------------------------------------------------------------------

		for p1 in ${phwh[@]}; do                             # -homozyg-window-het
			for p2 in ${phwm[@]}; do                         # -homozyg-winodow-missing
				for p3 in ${phws[@]}; do                     # -homozyg-window-snp
					for p4 in ${phzd[@]}; do                 # -homozyg-density
						for p5 in ${phzg[@]}; do             # -homozyg-gap
							for p6 in ${phwt[@]}; do         # -homozyg-window threshold
								for p7 in ${phzs[@]}; do     # -homozyg-snp
									for p8 in ${phzk[@]}; do # -homozyg-kb

										PARAM_SUFFIX=_phwh_${p1}_phwm_${p2}_phws_${p3}_phzd_${p4}_phzg_${p5}_phwt_${p6}_phzs_${p7}_phzk_${p8}
										IN_FILE=${d}_sample_cvg_${c}_plink
										OUT_FILE=${d}_sample_cvg_${c}_plink_roh${PARAM_SUFFIX}

										start_logging "PLINK ROH ID - ${PARAM_SUFFIX}"

										plink --bim ${OUTPUT_DIR}/${IN_FILE}.bim \
											--bed ${OUTPUT_DIR}/${IN_FILE}.bed \
											--fam ${OUTPUT_DIR}/${IN_FILE}.fam \
											--homozyg-window-het ${p1} \
											--homozyg-window-missing ${p2} \
											--homozyg-window-snp ${p3} \
											--homozyg-density ${p4} \
											--homozyg-gap ${p5} \
											--homozyg-window-threshold ${p6} \
											--homozyg-snp ${p7} \
											--homozyg-kb ${p8} \
											--out ${OUTPUT_DIR}/${OUT_FILE}

										stop_logging
									done
								done
							done
						done
					done
				done
			done
		done
	done
done