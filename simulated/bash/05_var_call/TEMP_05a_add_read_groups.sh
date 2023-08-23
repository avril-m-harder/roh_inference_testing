#!/bin/bash
#
#SBATCH --job-name=05a_run_downsample
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 03:30:00
#SBATCH --mem=20000
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=avrilharder@gmail.com
#SBATCH --array=0-4


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
# SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")
SCRIPT=05a_add_read_groups.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

declare -a dems=(decline large-1000)


# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load picard/2.23.9

# -----------------------------------------------------------------------------
# Run Haplotype caller on all sample files
# -----------------------------------------------------------------------------

for d in ${dems[@]}; do

	SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_ID_list_${d}.txt

	# Create output directory for each coverage level.

	CVG_OUTPUT_DIR=${OUTPUT_DIR}/${d}_sample_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}
	mkdir ${CVG_OUTPUT_DIR}

	# Set input directory for this coverage

	CVG_INPUT_DIR=${INPUT_DIR}/${d}_sample_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}

	# Process each individual sample file ----------------------------------

	while read -a line; do

		# # ------------------------------------------------------------------
		# # Add read group information
		# # ------------------------------------------------------------------

		IN_FILE=${d}_${line[0]}_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}.bam
		OUT_FILE=${d}_${line[0]}_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}_rgroups.bam
		RGROUPS_FILE=${OUT_FILE}
		echo IN_FILE: ${IN_FILE}
		echo OUT_FILE: ${OUT_FILE}
		start_logging "picard/AddOrReplaceReadGroups - ${OUT_FILE}"

		## Do picard read groups

		# Format for Alabama Supercomputer
		# java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar \

		# Format for Easley.auburn.edu
		java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
			I=${CVG_INPUT_DIR}/${IN_FILE} \
			O=${CVG_OUTPUT_DIR}/${OUT_FILE} SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=${line[0]} RGPU=unit4

		stop_logging

		# ------------------------------------------------------------------
		# Index BAM files
		# ------------------------------------------------------------------

		OUT_FILE=${d}_${line[0]}_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}_rgroups.bai
		echo BAM_OUT_FILE: ${OUT_FILE}
		start_logging "picard/BuildBamIndex - ${OUT_FILE}"

		## Run picard BuildBamIndex

		# Format for Alabama Supercomputer
		# java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/BuildBamIndex.jar \

		# Format for Easley.auburn.edu
		java -jar /tools/picard-2.23.9/libs/picard.jar BuildBamIndex \
			I=${CVG_OUTPUT_DIR}/${RGROUPS_FILE}

		stop_logging

	done <${SAMPLE_ID_LIST}

done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
