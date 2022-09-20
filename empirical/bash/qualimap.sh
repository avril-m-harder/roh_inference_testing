#!/bin/bash

#SBATCH --job-name=qualimap
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=80000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std (alternative = general)

#### Link to Qualimap command line overview:
## http://qualimap.conesalab.org/doc_html/command_line.html#multi-sample-bam-qc


##  Set username
USER=avrilh

## Set project name
PROJ=qualimap

#### IF aubaxh002_qualimap DIRECTORY HAS BEEN DELETED ####
## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/
# 
# ## Set permissions for directory
# chmod 700 /scratch/${USER}/${PROJ}/
# 
# ## cd into working scratch directory
# cd /scratch/${USER}/${PROJ}/
# 
# ## download and unzip qualimap
# wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip
# unzip qualimap_v2.2.1.zip
# cd qualimap_v2.2.1/


#### IF ALREADY DOWNLOADED ####
## cd into working scratch directory
cd /scratch/${USER}/${PROJ}/qualimap_v2.2.1/

## --------------------------------
## Load modules
module load R


## --------------------------------
## Run Qualimap in Multi-sample BAM QC mode, providing raw BAM files
# ./qualimap multi-bamqc \
# -d /home/amh0254/roh_param_project/sample_lists/qualimap_bam_files.txt \
# --run-bamqc \
# --java-mem-size=16G \
# --paint-chromosome-limits \
# -outformat PDF \
# -outfile multibamqc_report_01.pdf


## Run Qualimap on just a single sample for a quick methods check
# ./qualimap bamqc \
# -bam /scratch/avrilh/rohparam_03_gatk/ERR1474983_rgroups_dupmarked_fixmate.bam \
# --java-mem-size=16G \
# --paint-chromosome-limits \
# -outdir qualimap_out \
# -outformat PDF \
# -outfile multibamqc_report_01.pdf


## Run Qualimap on 15 tas devil samples
while read -a line
do

./qualimap bamqc \
-bam /scratch/avrilh/rohparam_02_read_mapping/sorted_bams/${line[0]}_rgroups.bam \
--java-mem-size=80G \
-nt 20 \
--paint-chromosome-limits \
-outdir qualimap_out \
-outformat PDF \
-outfile bamqc_report_${line[0]}.pdf

done < /home/amh0254/roh_param_project/sample_lists/15_samples.txt