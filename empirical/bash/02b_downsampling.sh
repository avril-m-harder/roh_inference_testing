#!/bin/bash

#SBATCH --job-name=downsamp
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=80000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std

##  Set username
# USER=avrilh

## Set project name
# PROJ=rohparam_03_gatk

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
# chmod 700 /scratch/${USER}/${PROJ}/

## cd into directory
cd /scratch/avrilh/rohparam_02_read_mapping/sorted_bams/


## --------------------------------
## Load modules
module load samtools/1.11


## --------------------------------
## For each of the 15 samples, downsample to roughly 5X,  10X, 15X, 30X, and 50X, 
## if necessary.

while read -a line
do

samtools view -s 0.11 -@ 19 \
	-o ${line[0]}_rgroups_5X.bam \
	${line[0]}_rgroups.bam
	
samtools view -s 0.23 -@ 19 \
	-o ${line[0]}_rgroups_10X.bam \
	${line[0]}_rgroups.bam
	
samtools view -s 0.34 -@ 19 \
	-o ${line[0]}_rgroups_15X.bam \
	${line[0]}_rgroups.bam
	
samtools view -s 0.68 -@ 19 \
	-o ${line[0]}_rgroups_30X.bam \
	${line[0]}_rgroups.bam
	
done < /home/amh0254/roh_param_project/sample_lists/15_samples.txt
