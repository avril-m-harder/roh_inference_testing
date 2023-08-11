#!/bin/bash
#
#SBATCH --job-name=_demo__02b_concat_fastq
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mem=4000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=02_read_sim_and_qc
PREV_STEP=01_slim
SCRIPT=02b_concat_fastq.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Concatenate forward and reverse read files to a single file for each
# individual
#
# OUTPUT_DIR and SAMPLE_ID_LIST are set in init_script_vars.sh
# -----------------------------------------------------------------------------

# cd to project directory
cd ${OUTPUT_DIR}

# Concatenate files

# for d in ${dems[@]}; do

	SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_ID_list__demo_.txt
	
	while read -a line; do

		start_logging "Concat fastq - ${line[0]}"

		cat _demo__${line[0]}_11.fq _demo__${line[0]}_21.fq > _demo__${line[0]}_f.fq
		cat _demo__${line[0]}_12.fq _demo__${line[0]}_22.fq > _demo__${line[0]}_r.fq

		stop_logging

	done < ${SAMPLE_ID_LIST}

# done
