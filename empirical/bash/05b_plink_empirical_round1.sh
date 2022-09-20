#!/bin/bash

#SBATCH --job-name=plink
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 300:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load plink/1.9

# -----------------------------------------------------------------------------
# Copy input VCF files to scratch
# -----------------------------------------------------------------------------

cp /scratch/avrilh/rohparam_03_gatk/sorted_tasdev_all_chroms_15samps_*.recode.vcf \
/scratch/avrilh/rohparam_05c_plink_empirical

cd /scratch/avrilh/rohparam_05c_plink_empirical

# ----------------------------------------------------------------------
# Convert VCF to PLINK format(s)
# ----------------------------------------------------------------------

declare -a cvgX=(fullcovg 5X 10X 15X 30X)

for c in ${cvgX[@]}; do
	
	IN_FILE=sorted_tasdev_all_chroms_15samps_${c}.recode.vcf
	OUT_FILE=tasdev_cvg_${c}_plink

	start_logging "Convert VCF to PLINK"

	plink \
		--vcf ${IN_FILE} \
		--allow-extra-chr \
		--out ${OUT_FILE}

	stop_logging

done

# -----------------------------------------------------------------------------
# PLINK parameters - define arrays for each of the PLINK parameters we want to
# vary and test.
# -----------------------------------------------------------------------------

declare -a phwh=(0 1 2)                   # Values for -homozyg-window-het
declare -a phwm=(2 5 50)                  # Values for -homozyg-window-missing
declare -a phws=(50 100 1000)             # Values for -homozyg-window-snp
declare -a phzd=(50)                      # Values for -homozyg-density
declare -a phzg=(500 1000)                # Values for -homozyg-gap
declare -a phwt=(0.01 0.05 0.1)           # Values for -homozyg-window-threshold
declare -a phzs=(10 100 1000)             # Values for -homozyg-snp
declare -a phzk=(100)                     # Values for -homozyg-kb

declare -a cvgX=(fullcovg 5X 10X 15X 30X)

# ----------------------------------------------------------------------
# Run PLINK to ID ROHs
# ----------------------------------------------------------------------

for c in ${cvgX[@]}; do
	for p1 in ${phwh[@]}; do                             # -homozyg-window-het
		for p2 in ${phwm[@]}; do                         # -homozyg-winodow-missing
			for p3 in ${phws[@]}; do                     # -homozyg-window-snp
				for p4 in ${phzd[@]}; do                 # -homozyg-density
					for p5 in ${phzg[@]}; do             # -homozyg-gap
						for p6 in ${phwt[@]}; do         # -homozyg-window threshold
							for p7 in ${phzs[@]}; do     # -homozyg-snp
								for p8 in ${phzk[@]}; do # -homozyg-kb

									PARAM_SUFFIX=_phwh_${p1}_phwm_${p2}_phws_${p3}_phzd_${p4}_phzg_${p5}_phwt_${p6}_phzs_${p7}_phzk_${p8}
									IN_FILE=tasdev_cvg_${c}_plink
									OUT_FILE=tasdev_cvg_${c}_plink_roh${PARAM_SUFFIX}

									start_logging "PLINK ROH ID - ${PARAM_SUFFIX}"

									plink --bim ${IN_FILE}.bim \
										--bed ${IN_FILE}.bed \
										--fam ${IN_FILE}.fam \
										--allow-extra-chr \
										--homozyg-window-het ${p1} \
										--homozyg-window-missing ${p2} \
										--homozyg-window-snp ${p3} \
										--homozyg-density ${p4} \
										--homozyg-gap ${p5} \
										--homozyg-window-threshold ${p6} \
										--homozyg-snp ${p7} \
										--homozyg-kb ${p8} \
										--out ${OUT_FILE}

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