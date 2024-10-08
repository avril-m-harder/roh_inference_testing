#!/bin/bash
#SBATCH --job-name=_demo__02c_read_sim_and_qc
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 02:00:00
#SBATCH --mem=4000
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=02_read_sim_and_qc
PREV_STEP=01_slim
SCRIPT=02c_fastqc.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh


## --------------------------------
## Load modules# module load trimmomatic/0.38
module load fastqc/0.11.9

## Run fastqc on 2 read files per sample

# for d in ${dems[@]}; do

	# Create directory to save fastqc reports
	cd ${OUTPUT_DIR}
	cd ..
	FASTQC_OUT_DIR=fastqc_reports__demo_
	mkdir ${FASTQC_OUT_DIR}

	SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_ID_list__demo_.txt

	while read -a line; do

		start_logging "fastqc - ${line[0]}"

		fastqc -t 5 -o ./${FASTQC_OUT_DIR} \
			${OUTPUT_DIR}/_demo__${line[0]}_f.fq \
			${OUTPUT_DIR}/_demo__${line[0]}_r.fq

		stop_logging

	done < ${SAMPLE_ID_LIST}

# done

cd /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/02_read_sim_and_qc
mv slurm*.out ./script_stdouts/