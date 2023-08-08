#!/bin/bash

#SBATCH --job-name=01_run_slim
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 300:00:00
#SBATCH --mem=12000
#SBATCH --output=01_slim-%j.out
#SBATCH --error=error-01_slim-%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=avrilharder@gmail.com


# IMPORTANT: before running this script, you must set SLiM parameters in
# init_script_vars.sh


# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=01_slim
SCRIPT=01_run_slim

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load bcftools/1.17
module load samtools/1.17
module load bedtools/2.30.0

# -----------------------------------------------------------------------------
# [NOTE: SLiM is now run on local machine (see README in 
# simulated/stoffel_et_al_sheep_roh_AMH/), and resulting VCF files are uploaded
# to cluster for further processing below ]
#
# Prep ancestral fasta file from tas dev reference chromosome 1. start and end
# coordinates generated randomly in R (local)
# -----------------------------------------------------------------------------

mkdir ${OUTPUT_DIR}/${SLIM_OUT_DIR}
cd ${OUTPUT_DIR}/${SLIM_OUT_DIR}
cp ${HOME_STEP_DIR}/*.vcf .

## issues with bedtools not recognizing chrom name, just rename it.
echo -e "1\t129649026\t159649025\t1\n" \
> subset.bed

cp /scratch/avrilh/archive_first_roh_param_runs/rohparam_01b_assemb/tasdev_genbank_assem.fna .
sed -i 's/LR735554.1 Sarcophilus harrisii genome assembly, chromosome: 1/1/g' \
tasdev_genbank_assem.fna

bedtools getfasta \
-fo ${REF_GENOME_FILE_PATH}/${REF_GENOME_FILE_NAME}.fasta \
-nameOnly \
-fi tasdev_genbank_assem.fna \
-bed subset.bed

rm tasdev_genbank_assem.fna*


# -----------------------------------------------------------------------------
# Format SLiM output for read simulation - FILE_LABELS and SAMPLE_ID_LIST and
# FASTA_OUT_DIR are set in init_script_vars.sh. FASTA_OUT_DIR is used as the
# input directory in 02a_run_art.sh
# -----------------------------------------------------------------------------

start_logging "Format SLiM Output for read simulation - ${SLIM_OUT_DIR}"

VCF_OUT_DIR=${OUTPUT_DIR}/sample_vcf_files_${FILE_LABELS}
mkdir ${VCF_OUT_DIR}
VCF_FILE_LIST=${OUTPUT_DIR}/vcf_file_list_${FILE_LABELS}.txt

## Compress output VCF
bgzip ${OUTPUT_DIR}/${SLIM_OUT_DIR}/final_pop.vcf

## Index compressed VCF
tabix -f ${OUTPUT_DIR}/${SLIM_OUT_DIR}/final_pop.vcf.gz

## Split output VCF into sample-specific VCF files

bcftools +split -O z -o ${VCF_OUT_DIR}/ ${OUTPUT_DIR}/${SLIM_OUT_DIR}/final_pop.vcf.gz

## Randomly select sample VCF files for conversation to FASTAs
## >> Selecting 100, can be downsampled later

ls ${VCF_OUT_DIR} i*.vcf.gz | sort -R | tail -100 >${VCF_FILE_LIST}
sed 's/.vcf.gz//g' ${VCF_FILE_LIST} >${SAMPLE_ID_LIST}

## Convert sample VCF files to two separate haplotype FASTAs per individual

while read -a line; do
#     tabix ${VCF_OUT_DIR}/${line[0]}.vcf.gz

    bcftools norm --check-ref s \
        --fasta-ref ${REF_GENOME_FILE} \
        --multiallelics - \
        --do-not-normalize \
        --output ${VCF_OUT_DIR}/norm_${line[0]}.vcf \
        /${VCF_OUT_DIR}/${line[0]}.vcf.gz

    bgzip ${VCF_OUT_DIR}/norm_${line[0]}.vcf
    tabix ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

    bcftools +fixref ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz -- \
        --fasta-ref ${REF_GENOME_FILE}

    bcftools consensus --haplotype 1 \
        --fasta-ref ${REF_GENOME_FILE} \
        --output ${FASTA_OUT_DIR}/${line[0]}_1.fasta \
        ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

    bcftools consensus --haplotype 2 \
        --fasta-ref ${REF_GENOME_FILE} \
        --output ${FASTA_OUT_DIR}/${line[0]}_2.fasta \
        ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

done <${SAMPLE_ID_LIST}

stop_logging

mail -s 'SLiM run finished - submit Rviz' ${EMAIL} <<<'SLiM run finished'

mv *.out /home/amh0254/roh_param_project/roh_inference_testing/simulated/script_stdouts
mv *.err /home/amh0254/roh_param_project/roh_inference_testing/simulated/script_stdouts

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/99_includes/backup_output.sh
