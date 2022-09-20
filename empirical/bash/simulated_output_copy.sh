#!/bin/bash

#SBATCH --job-name=output_copy
#SBATCH --partition=jrw0107_std 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

cp /scratch/kbk0024/roh_param_project/data/06_roh_id/output/*.hom \
/home/amh0254/roh_param_project/simulated_output/plink/
