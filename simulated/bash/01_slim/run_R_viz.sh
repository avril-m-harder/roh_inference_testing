#!/bin/bash
#SBATCH --job-name=run_r_viz
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mem=8000
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=avrilharder@gmail.com
#SBATCH --output=01_Rviz-%j.out
#SBATCH --error=error-01_Rviz-%j.err


## --------------------------------
## Load R + run script 
module load R/4.2.3

Rscript \
/home/amh0254/roh_param_project/roh_inference_testing/simulated/bash/01_slim/\
slim_output_roh_viz_Easley.R

# mail -s 'R viz finished' avrilharder@gmail.com <<< 'R viz finished'

mv *.out /home/amh0254/roh_param_project/roh_inference_testing/simulated/script_stdouts
mv *.err /home/amh0254/roh_param_project/roh_inference_testing/simulated/script_stdouts
