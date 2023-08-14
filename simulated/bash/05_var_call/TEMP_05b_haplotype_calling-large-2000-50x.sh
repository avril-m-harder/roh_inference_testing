#!/bin/bash

#SBATCH --job-name=05b_haplotype_calling
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH --mem=8000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com


#  +----------------------+
#  | REQUEST 4 CPU + 60gb |
#  +----------------------+

#  Before running this script, make sure that dictionary and index files exist
#  for the reference genome file. To create dictionary file:
#
#     gatk-launch CreateSequenceDictionary -R ref.fasta
#
#  to create index file:
#
#     samtools faidx ref.fasta
#
#
#  NOTES for running on Alabama Supercomputer:
#
#  NOTE: After all of the scripts finish running, check the output directories -
#        sample_cvg_XX and make sure they have the expected number of files.
#        Sometimes one or two samples don't get processed. If there are fewer
#        files than expected, compare the list of files with the list of
#        samples in SAMPLE_ID_LIST to determine which samples are missing. Then
#        re-run the matching 05d_gtyp_pXX_c_YYx.sh script in 05_var_call
#        scripts directory.

#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05b_haplotype_calling.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load gatk/4.1.9.0

# -----------------------------------------------------------------------------
# Run Haplotype caller on all sample files
# -----------------------------------------------------------------------------


SAMPLE_ID_LIST=${INIT_OUTPUT_DIR}/sample_ID_list_large-2000.txt

# Create output directory for each coverage level.

CVG_OUTPUT_DIR=${OUTPUT_DIR}/large-2000_sample_cvg_50x

# Process each individual sample file ----------------------------------

while read -a line; do

	# ------------------------------------------------------------------
	# Run HaplotypeCaller in GVCF mode
	# ------------------------------------------------------------------

	OUT_FILE=large-2000_${line[0]}_cvg_50x.g.vcf
	RGROUPS_FILE=large-2000_${line[0]}_cvg_50x_rgroups.bam
	start_logging "gatk Haplotype Caller - ${OUT_FILE}"

	# Run gatk HaplotypeCaller

	gatk HaplotypeCaller \
		-R ${REF_GENOME_FILE} \
		-I ${CVG_OUTPUT_DIR}/${RGROUPS_FILE} \
		-stand-call-conf 20.0 \
		--emit-ref-confidence GVCF \
		-O ${CVG_OUTPUT_DIR}/${OUT_FILE}

	stop_logging

done <${SAMPLE_ID_LIST}


# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
