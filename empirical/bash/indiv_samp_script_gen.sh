#!/bin/bash

SCRIPT=/home/amh0254/roh_param_project/scripts/03a_preBQSR_snp_calling.sh
SCRATCHDIR=/scratch/avrilh/rohparam_03_gatk

while read -a line
do

#sed "s/_samp_/${line[0]}/g" ${SCRIPT} > ${SCRATCHDIR}/03a_preBQSR_snp_calling_${line[0]}.sh

## check scripts before submitting
sbatch ${SCRATCHDIR}/03a_preBQSR_snp_calling_${line[0]}.sh

done < /home/amh0254/roh_param_project/sample_lists/15_samples.txt
