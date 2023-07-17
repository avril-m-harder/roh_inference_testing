#!/bin/bash

#SBATCH --job-name=02_read_sim_and_qc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 07:00:00
#SBATCH --mem=4000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com

#  This script must be made executable like this
#    chmod +x my_script
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=02_read_sim_and_qc
PREV_STEP=01b_create_subsample_lists
SCRIPT=02a_run_art.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load artmountrainier/2016.06.05

# -----------------------------------------------------------------------------
# Run read simulation
#
# INPUT_DIR, SLIM_OUT_DIR, FASTA_OUT_DIR, OUTPUT_DIR, and SAMPLE_ID_LIST are
# set in  init_script_vars.sh
# -----------------------------------------------------------------------------

# Starting with 50X coverage for now (25X per haploid genome, i.e., FASTA file)
# Will downsample unfiltered BAM files generated by BWA.
# 150-bp paired-end reads (--len and --paired parameters).
# Using the HiSeq 2500 error model (--seqSys parameter).

while read -a line; do

    for a in 1 2; do

        start_logging "ART Read - ${line[0]}_$a.fasta"

        art_illumina \
            --seqSys HS25 \
            -i ${FASTA_OUT_DIR}/${line[0]}_$a.fas \
            --paired \
            -na \
            --len 150 \
            --fcov 25 \
            -m 500 \
            -s 75 \
            -o ${OUTPUT_DIR}/${line[0]}_$a

        # echo ${FASTA_OUT_DIR}/${line[0]}_$a.fas
        # echo ${OUTPUT_DIR}/${line[0]}_$a

        stop_logging

    done

done <${SAMPLE_ID_LIST}

# mail -s 'ART run finished' avrilharder@gmail.com <<<'SLiM run finished'

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/backup_output.sh
