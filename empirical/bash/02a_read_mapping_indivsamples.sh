#!/bin/bash

#SBATCH --job-name=02_roh__samp_
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's partition: jrw0107_std

##  Set username
USER=avrilh

## Set project name
PROJ=rohparam_02_read_mapping

## Set reference genome
REF=/scratch/avrilh/rohparam_01b_assemb/tasdev_genbank_assem.fna.gz

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


## --------------------------------
## Load modules
module load bwa/0.7.17
module load samtools/1.11
module load picard/2.23.9
#module load bedtools/2.29.2


## --------------------------------
## Align reads to genome
# while read -a line
# do

bwa mem -t 20 -M ${REF} \
/scratch/avrilh/rohparam_01a_dl_and_qc/_samp__1_val_1.fq.gz \
/scratch/avrilh/rohparam_01a_dl_and_qc/_samp__2_val_2.fq.gz \
> _samp_.bam


## add read group information 
java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
I=_samp_.bam \
O=_samp__rgroups.bam \
SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=_samp_ RGPU=unit1

## remove original BAM file or scratch gets completely filled up --
## check read #s against SRA metadata
rm _samp_.bam

## get some mapping stats
samtools flagstat -@ 20 -O tsv _samp__rgroups.bam > _samp__rgroups_flagstat_out.txt


## --------------------------------
cp _samp__rgroups_flagstat_out.txt /home/amh0254/roh_param_project/02_read_mapping/flagstat_out/

