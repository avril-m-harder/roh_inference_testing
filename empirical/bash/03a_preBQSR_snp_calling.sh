#!/bin/bash

#SBATCH --job-name=roh__samp__03a_gatk
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=80000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std


##  Set username
USER=avrilh

## Set project name
PROJ=rohparam_03_gatk

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into directory
cd /scratch/${USER}/${PROJ}/

ref=/scratch/avrilh/rohparam_01b_assemb/tasdev_genbank_assem.fna

## --------------------------------
## Load modules
module load samtools/1.11
module load picard/2.23.9
module load gatk/4.1.9.0


## --------------------------------
## Mark duplicate reads - GATK tools ignore them, no need to remove, just flag
for c in 5 10 15 30
do

java -jar /tools/picard-2.23.9/libs/picard.jar MarkDuplicates \
I=/scratch/avrilh/rohparam_02_read_mapping/sorted_bams/_samp__rgroups_${c}X.bam \
O=_samp__rgroups_${c}X_dupmarked.bam \
M=_samp__rgroups_${c}X_markdup_metrics.txt \
MAX_RECORDS_IN_RAM=250000

## Fix mate information
java -jar /tools/picard-2.23.9/libs/picard.jar FixMateInformation \
I=_samp__rgroups_${c}X_dupmarked.bam \
O=_samp__rgroups_${c}X_dupmarked_fixmate.bam

## Index BAM files -- .bai doesn't work for the 3 largest chroms (need .csi indexing)
# java -jar /tools/picard-2.23.9/libs/picard.jar BuildBamIndex \
# I=_samp__rgroups_dupmarked_fixmate.bam

samtools index -c -@ 3 _samp__rgroups_${c}X_dupmarked_fixmate.bam

## Run HaplotypeCaller in GVCF mode
gatk HaplotypeCaller \
	-R $ref \
	-I _samp__rgroups_${c}X_dupmarked_fixmate.bam \
	-stand-call-conf 20.0 \
	--emit-ref-confidence GVCF \
	-O _samp__${c}X_nobaseQrecal.g.vcf
	
done