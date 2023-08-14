#!/bin/bash

#SBATCH --job-name=_demo__03_read_align
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 20:00:00
#SBATCH --mem=60000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=03_read_align
PREV_STEP=02_read_sim_and_qc
SCRIPT=03_run_bwa.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load bwa/2.0
module load samtools/1.17

# -----------------------------------------------------------------------------
# Index reference genomes
# -----------------------------------------------------------------------------

# TO DO: move this indexing to step O1  - make it step 01c_index_ref_genome
# For now, just ran it manually since using makefile

# bwa index ${REF_GENOME_FILE}

# -----------------------------------------------------------------------------
# Align reads to reference genome
# -----------------------------------------------------------------------------

# SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_ID_list__demo_.txt

cd ${OUTPUT_DIR}
# mkdir _demo_; cd _demo_
# 
# while read -a line; do
# 
# 	OUT_FILE=_demo__${line[0]}_genome.bam
# 	start_logging "bwa align - ${OUT_FILE}"
# 
# 	bwa mem -t 20 -M \
# 		${REF_GENOME_FILE} \
# 		${INPUT_DIR}/_demo_/_demo__${line[0]}_f.fq \
# 		${INPUT_DIR}/_demo_/_demo__${line[0]}_r.fq \
# 		> ${OUT_FILE}
# 
# 	stop_logging
# 
# done < ${SAMPLE_ID_LIST}

# -------------------------------------------------------------------------
# Sort aligned read bam files
# -------------------------------------------------------------------------

# while read -a line; do
# 
# 	OUT_FILE=_demo__${line[0]}_genome_sorted.bam
# 	start_logging "samtools sort - ${OUT_FILE}"
# 
# 	samtools sort -@ 19 \
# 		-o ${OUT_FILE} \
# 		_demo__${line[0]}_genome.bam
# 
# 	stop_logging
# 
# done < ${SAMPLE_ID_LIST}
	

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

cp ./_demo_/_demo__*_genome_sorted.bam \
/home/amh0254/roh_param_project/roh_inference_testing/simulated/data/03_read_align