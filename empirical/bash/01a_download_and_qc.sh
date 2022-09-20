#!/bin/bash

#SBATCH --job-name=01a_dl_qc__group_
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com
#SBATCH --output=/home/amh0254/roh_param_project/scripts/stdouts/%j.out
#SBATCH --error=/home/amh0254/roh_param_project/scripts/stdouts/error-%j.err

##  Set username
USER=avrilh

## Set project name
PROJ=rohparam_01a_dl_and_qc

## Set batch group
GROUP=_group_

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/_group_


## --------------------------------
## Activate conda environment and load other modules 
cd
source activate env1
module load python/3.8.6
module load fastqc/0.11.9
module load sratoolkit/2.11.0

## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/_group_/


## --------------------------------
## Download FASTQ files from SRA
# cat /home/amh0254/roh_param_project/sample_lists/_group_.txt | while read line
# 	do
# 	fastq-dump --split-files --origfmt --gzip ${line}
# 	done

## --------------------------------
## Run fastqc on all 4 read files per sample
# mkdir pretrim_output/
# mkdir posttrim_output/
# mkdir stdouts/

# while read -a line
# do
# 	fastqc -t 1 -o pretrim_output/ \
# 	${line[0]}_1.fastq.gz \
# 	${line[0]}_2.fastq.gz
# done < /home/amh0254/roh_param_project/sample_lists/_group_.txt
# 
# cp pretrim_output/* /home/amh0254/roh_param_project/01_read_qc_trimming/pretrim_fastqc/

## Run TrimGalore to clean and trim adapters from reads,
## then run FastQC on the trimmed files
while read -a line
do

	/home/amh0254/programs/trimgalore/TrimGalore-0.6.6/trim_galore \
	--paired --quality 20 --phred33 \
	--fastqc_args "--outdir posttrim_output/" \
	--length 30 \
	${line[0]}_1.fastq.gz \
	${line[0]}_2.fastq.gz

done < /home/amh0254/roh_param_project/sample_lists/_group_.txt

cp posttrim_output/* /home/amh0254/roh_param_project/01_read_qc_trimming/posttrim_fastqc/

## Deactivate conda environment
conda deactivate 

## Clean up
# mv *.out stdouts/
