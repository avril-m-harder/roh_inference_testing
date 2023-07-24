#!/bin/bash
#SBATCH --job-name=run_r_viz
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mem=8000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=avrilharder@gmail.com


## --------------------------------
## Load R + run script 
module load R/4.2.3

Rscript \
/home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/01_slim/\
slim_output_roh_viz_Easley.R

mail -s 'R viz finished' avrilharder@gmail.com <<< 'R viz finished'
