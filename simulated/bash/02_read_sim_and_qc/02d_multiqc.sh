#!/bin/bash
#SBATCH --job-name=02d_multiqc
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mem=4000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=02_read_sim_and_qc
PREV_STEP=01_slim
SCRIPT=02d_multiqc.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh


## --------------------------------
## Load modules
module load python/3.9.2

# pip install --user multiqc

## Run multiQC in each demo scenario's fastQC folder

for d in ${dems[@]}; do

	FASTQC_OUT_DIR=fastqc_reports_${d}
	cd ${OUTPUT_DIR}
	cd ../${FASTQC_OUT_DIR} 
	multiqc .
	
	cp multiqc_report.html ${HOME_STEP_DIR}
	mv ${HOME_STEP_DIR}/slurm*.out ${HOME_STEP_DIR}/script_stdouts/

done
