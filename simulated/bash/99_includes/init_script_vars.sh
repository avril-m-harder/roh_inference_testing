#!/bin/bash
#
# initialize_script_variables.sh
#
# This script initializes variables that are used in multiple scripts in the
# ROH calling workflow. It should be included in all scripts for each step of
# the workflow. Put the following line in the scripts:
#
#     source /scratch/aubkbk001/roh_param_project/init_script_vars.sh
#
# replacing aubkbk001 with your username.
#
# 2022-04-28 - Modify log file creation code - see if this script is being
#              run as an array job. If so, add the job id and task id to the
#              log file name

##  Set username
USER1=amh0254
USER2=avrilh

## Set project name
PROJECT=roh_param_project/roh_inference_testing

## Set base directories
USER_HOME_DIR=/home/${USER1}
BASE_DIR=/scratch/${USER2}

## Set email address
EMAIL=avrilharder@gmail.com

## Set initial step name
INIT_STEP=01_slim

# -----------------------------------------------------------------------------
# Create working directories
# -----------------------------------------------------------------------------

## Create variable pointing to the base directory for this project and step

HOME_STEP_DIR=${USER_HOME_DIR}/${PROJECT}/simulated/data/${STEP}/
BASE_STEP_DIR=${BASE_DIR}/${PROJECT}/data/${STEP}/

## Set the input directory, which is the output directory from the previous
## step
INPUT_DIR=${BASE_DIR}/${PROJECT}/simulated/data/${PREV_STEP}/output

## Create a working directory on /scratch
WORK_DIR=${BASE_DIR}/${PROJECT}/simulated/data/${STEP}
mkdir -p ${WORK_DIR}
chmod -R 700 ${WORK_DIR}

## Create output directory in working directory
OUTPUT_DIR=${WORK_DIR}/output
mkdir ${OUTPUT_DIR}
chmod -R 700 ${OUTPUT_DIR}

# -----------------------------------------------------------------------------
## Set Slim parameters
# -----------------------------------------------------------------------------

# directory and file paths

SLIM_OUT_DIR=slim_output

INIT_OUTPUT_DIR=${BASE_DIR}/${PROJECT}/simulated/data/${INIT_STEP}/output/

REF_GENOME_FILE_PATH=${INIT_OUTPUT_DIR}

REF_GENOME_FILE_NAME=ancestral

REF_GENOME_FILE=${REF_GENOME_FILE_PATH}/${REF_GENOME_FILE_NAME}.fasta

# BASE_SAMPLE_ID_LIST contails a list of all the samples we generated in SLiM.
# We won't necessarily use all of those elements.

# SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_id_list_${FILE_LABELS}.txt
# BASE_SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_id_list.txt
# 
# FASTA_OUT_DIR=${INIT_OUTPUT_DIR}/sample_fasta_files_${FILE_LABELS}
# FASTA_OUT_DIR=${INIT_OUTPUT_DIR}/slim_output_fasta_files
# mkdir ${FASTA_OUT_DIR}

# DEMO_SCEN_NAMES=${HOME_STEP_DIR}/demo_scenario_names.txt

# -----------------------------------------------------------------------------
# Create arrays defining coverage levels and population sizes. These arrays are
# used in:
#
#    04_downsample.sh
#    05a_snp_calling.sh
#    05b_genomicsdb.sh
#    06a_run_bcftoolsROH.sh
#    06b_run_plink.sh
# -----------------------------------------------------------------------------

## Create arrays of downsample levels to be used.
## Coverage level in NNx for display in output file names
declare -a cvgX=(50x 30x 15x 10x 05x)

## Coverage level fraction to supply to samtools
declare -a cvgP=(1.0 0.6 0.3 0.2 0.1)

## Names of demographic scenarios to be analyzed
declare -a dems=(small decline bottle large-2000)

## Get length of the coverage level arrays. Subtract 1 because arrays are zero
## based, and we'll iterate over the arrays from 0 to cvgCnt
cvgCnt=${#cvgX[@]}
let cvgCnt-=1

## Create array of population sizes we want to test. Populations MUST be listed
## in order of decreasing size. 01b_create_subsample_lists.sh depends on the
## populations being in descending order.

declare -a popN=(100 50 30)

popCnt=${#popN[@]}
let popCnt-=1

# SAMPLE_ID_LIST, set below, contains the list of the elements of the largest
# subsample we extract from the base list. e.g. if we create subsamples of 100,
# 50, and 30 individuals, SAMPLE_ID_LIST will point to the list of 100
# individuals.

SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_id_list_pop_${popN[0]}.txt

# -----------------------------------------------------------------------------
# Create log file
# -----------------------------------------------------------------------------

# Set name of log file to track execution times

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")
echo job:${SLURM_ARRAY_JOB_ID}
echo task:${SLURM_ARRAY_TASK_ID}

# If the job is being run as a slurm array job, add job and task ids to the
# file name.
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    LOG_FILE=${OUTPUT_DIR}/${SCRIPT}_${TIMESTAMP}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_log.txt
else
    LOG_FILE=${OUTPUT_DIR}/${SCRIPT}_${TIMESTAMP}_log.txt
fi

# # If log file exists, save a copy.

# if [ -f "$LOG_FILE" ]; then
#     TIMESTAMP=$(date "+%Y%m%d-%H%M%S")
#     mv ${LOG_FILE} ${OUTPUT_DIR}/${SCRIPT}_log_${TIMESTAMP}.txt
# fi

# Write header to log file

printf "%-80s   %8s   %8s   %8s\n" "Action - Output" "Start" "End" "Duration" >${LOG_FILE}

# -----------------------------------------------------------------------------
# Define Functions
# -----------------------------------------------------------------------------

# Log file functions - used in all scripts

start_logging() {

    ACTION=$1
    START_TIME_HHMM=$(date +" %T")
    START_TIME=$(date +%s)
}

stop_logging() {

    END_TIME_HHMM=$(date +"%T")
    END_TIME=$(date +%s)
    ELAPSED=$(expr $END_TIME - $START_TIME)
    DURATION=$(date -u -d @${ELAPSED} +"%T")
    printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}
}
