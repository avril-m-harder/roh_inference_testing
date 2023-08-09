#!/bin/bash

#SBATCH --job-name=01_run_slim
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 300:00:00
#SBATCH --mem=64000
#SBATCH --output=01_slim-%j.out
#SBATCH --error=error-01_slim-%j.err


# -----------------------------------------------------------------------------
# Load modules and install Python packages
# -----------------------------------------------------------------------------

cd

module load python/3.9.2
module load slim/4.0.1

# pip install --user msprime
# pip install --user pyslim


# -----------------------------------------------------------------------------
# Copy dirs to scratch
# -----------------------------------------------------------------------------

cd /scratch/avrilh/roh_param_project/roh_inference_testing/simulated/
cp -r /home/amh0254/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH .
cd stoffel_et_al_sheep_roh_AMH
mkdir output


# -----------------------------------------------------------------------------
# Run SLiM simulations and quick viz check
# -----------------------------------------------------------------------------

module load R/4.2.3

if [ -f "output/demo_results_overview.pdf" ]; then
	rm output/demo_results_overview.pdf
fi

while read -a line; do

	mkdir output/${line[0]}
	cd output/${line[0]}
	mkdir muts; mkdir out; mkdir roh; mkdir slim_code; mkdir trees; mkdir vcfs
	
	Rscript ../../scripts/EASLEY_slim_sims_pipeline.R \
	${line[0]} ${line[1]} ${line[2]} ${line[3]} \
	${line[4]} ${line[5]} ${line[6]} ${line[7]} ${line[8]}
	
	cd ../../

done < demo_scenarios.txt

cd ../../scripts/

Rscript EASLEY_demo_results_quick_check.R
