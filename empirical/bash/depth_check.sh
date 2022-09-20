#!/bin/bash

#SBATCH --job-name=_id_
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

module load samtools/1.11

cd /scratch/avrilh/rohparam_02_read_mapping/sorted_bams/

# while read -a line
# do

# samtools coverage --output covg__id_.txt _id__rgroups.bam

# samtools depth -o covg__id_.txt _id__rgroups.bam

# done < /home/amh0254/roh_param_project/sample_lists/15_samples.txt

for c in 5 10 15 30
	do
	samtools coverage -o ${c}X__id_.txt _id__rgroups_${c}X.bam
	done