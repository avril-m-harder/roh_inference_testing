#!/bin/bash

#SBATCH --job-name=01b_rohparam_assemb
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mem=8000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

##  Set username
USER=avrilh

## Set project name
PROJ=rohparam_01b_assemb

## Create directory for project
mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/


## --------------------------------
## Load modules
module load samtools/1.11
module load picard/2.23.9
module load bwa/0.7.17


## --------------------------------
## Download assembly as published
# wget \
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/635/505/GCA_902635505.1_mSarHar1.11/\
# GCA_902635505.1_mSarHar1.11_genomic.fna.gz
# 
# mv GCA_902635505.1_mSarHar1.11_genomic.fna.gz tasdev_genbank_assem.fna.gz


## --------------------------------
## Create indexes needed through GATK Best Practices
# bwa index tasdev_genbank_assem.fna.gz
# 
# samtools faidx tasdev_genbank_assem.fna.gz
# 
# java -jar /tools/picard-2.23.9/libs/picard.jar CreateSequenceDictionary \
# REFERENCE=tasdev_genbank_assem.fna.gz \
# OUTPUT=tasdev_genbank_assem.dictionary.bam

## --------------------------------
## Unzip FASTA and create indexes needed through GATK Best Practices
gunzip tasdev_genbank_assem.fna.gz

bwa index tasdev_genbank_assem.fna

samtools faidx tasdev_genbank_assem.fna

java -jar /tools/picard-2.23.9/libs/picard.jar CreateSequenceDictionary \
REFERENCE=tasdev_genbank_assem.fna \
OUTPUT=tasdev_genbank_assem.fna.dict
