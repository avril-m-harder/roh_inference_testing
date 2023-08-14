#!/bin/bash
#
#SBATCH --job-name=bottle_04_downsample
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH --mem=24000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=04_downsample
PREV_STEP=03_read_align
# This doesn't work when running via sbatch on easley. It does work when
# running the script in the shell manually, or running it under srun. Go figure.
# SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")
SCRIPT=04_run_downsample.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load samtools/1.17

# -----------------------------------------------------------------------------
# Do the downsampling.
# -----------------------------------------------------------------------------

# for d in ${dems[@]}; do

	SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_ID_list_bottle.txt

	for i in $(seq 0 $cvgCnt); do

		# # Create output directory for each coverage level.

		CVG_OUTPUT_DIR=${OUTPUT_DIR}/bottle_sample_cvg_${cvgX[i]}
		mkdir ${CVG_OUTPUT_DIR}

		# # Downsample each individual

		while read -a line; do

			OUT_FILE=bottle_${line[0]}_cvg_${cvgX[i]}.bam
			start_logging "samtools view - ${OUT_FILE}"

			samtools view -s ${cvgP[i]} -@ 19 \
				-o ${CVG_OUTPUT_DIR}/${OUT_FILE} \
				${INPUT_DIR}/bottle/bottle_${line[0]}_genome_sorted.bam
			
			samtools coverage ${CVG_OUTPUT_DIR}/${OUT_FILE}

			stop_logging

		done < ${SAMPLE_ID_LIST}

	done
	
# done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

