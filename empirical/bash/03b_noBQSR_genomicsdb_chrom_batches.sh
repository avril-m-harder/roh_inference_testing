#!/bin/bash

#SBATCH --job-name=roh_03b_gatk
#SBATCH --partition=jrw0107_std
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=80000
#SBATCH -t 300:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=avrilharder@gmail.com

## J's queue: jrw0107_std


##  Set username
USER=avrilh

## Set project name
PROJ=rohparam_03_gatk

## Create a directory on /scratch
# mkdir /scratch/${USER}/${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}/${PROJ}/

## cd into directory
cd /scratch/${USER}/${PROJ}/

ref=/scratch/avrilh/rohparam_01b_assemb/tasdev_genbank_assem.fna

## --------------------------------
## Load modules
module load samtools/1.11
module load picard/2.23.9
module load gatk/4.1.9.0
module load bcftools/1.11
module load vcftools/0.1.17


## --------------------------------
## Genotype GVCFs across all samples simultaneously --
## Only analyzes variants on autosomal chromosomes
for c in 5
do

	while read -a line
	do

	gatk --java-options "-Xmx80G" GenomicsDBImport \
		-V SRR9329973_$cX_nobaseQrecal.g.vcf \
		-V SRR9329974_$cX_nobaseQrecal.g.vcf \
		-V SRR9329975_$cX_nobaseQrecal.g.vcf \
		-V SRR9329976_$cX_nobaseQrecal.g.vcf \
		-V SRR9329977_$cX_nobaseQrecal.g.vcf \
		-V SRR9329978_$cX_nobaseQrecal.g.vcf \
		-V SRR9329979_$cX_nobaseQrecal.g.vcf \
		-V SRR9329980_$cX_nobaseQrecal.g.vcf \
		-V SRR9329981_$cX_nobaseQrecal.g.vcf \
		-V SRR9329982_$cX_nobaseQrecal.g.vcf \
		-V SRR9329983_$cX_nobaseQrecal.g.vcf \
		-V SRR9329986_$cX_nobaseQrecal.g.vcf \
		-V SRR9329989_$cX_nobaseQrecal.g.vcf \
		-V SRR9329990_$cX_nobaseQrecal.g.vcf \
		-V SRR9329993_$cX_nobaseQrecal.g.vcf \
		--genomicsdb-workspace-path ${line[0]}_$cX_database_gatk_genomics \
		--intervals ${line[0]}

	gatk --java-options "-Xmx80G" GenotypeGVCFs \
		-R $ref \
		-V gendb://${line[0]}_$cX_database_gatk_genomics \
		-O ${line[0]}_$cX_genotype_output_nobaseQrecal.vcf
	
	gatk --java-options "-Xmx80G" VariantFiltration \
		-R $ref \
		-V ${line[0]}_$cX_genotype_output_nobaseQrecal.vcf \
		-O ${line[0]}_$cX_genotype_output_nobaseQrecal_filtered.vcf \
		--filter-name "QD" \
		--filter-expression "QD < 2.0" \
		--filter-name "FS" \
		--filter-expression "FS > 60.0" \
		--filter-name "MQ" \
		--filter-expression "MQ < 40.0" \
		--filter-name "SOR" \
		--filter-expression "SOR > 5.0" \
		--filter-name "MQRankSum" \
		--filter-expression " MQRankSum < -3.0 || MQRankSum > 3.0" \
		--filter-name "ReadPosRankSum" \
		--filter-expression "ReadPosRankSum < -8.0"

	gatk --java-options "-Xmx80G" SelectVariants \
		-R $ref \
		-V ${line[0]}_$cX_genotype_output_nobaseQrecal_filtered.vcf \
		--select-type-to-include SNP \
		-select 'vc.isNotFiltered()' \
		-O filtered_${line[0]}_$cX_nobaseQrecal_SNPs.vcf


	## --------------------------------
	## Finally, remove SNPs within 5 bp of indels using bcftools before retaining only 
	## biallelic SNPs with SelectVariants

	bcftools filter --SnpGap 5 \
	--threads 4 \
	-o ${line[0]}_$cX_genotype_output_nobaseQrecal_filtered_nearindelfilt.vcf \
	${line[0]}_$cX_genotype_output_nobaseQrecal_filtered.vcf

	gatk --java-options "-Xmx80G" SelectVariants \
		-R $ref \
		-V ${line[0]}_$cX_genotype_output_nobaseQrecal_filtered_nearindelfilt.vcf \
		--select-type-to-include SNP \
		-select 'vc.isNotFiltered()' \
		-restrict-alleles-to BIALLELIC \
		-O final_filtered_${line[0]}_$cX_nobaseQrecal_SNPs_nearindelfilt.vcf
	
	echo ${line[0]} >> done_contigs.txt

	done < /home/amh0254/roh_param_project/sample_lists/tasdev_chroms.txt


	## --------------------------------
	## Concatenate files for all chromosomes into a single file
	ls final_filtered_*_$cX_nobaseQrecal_SNPs_nearindelfilt.vcf > \
	top_15_samps_autosomal_chroms.txt

	bcftools concat \
	-f top_15_samps_autosomal_chroms.txt \
	--threads 4 \
	-o final_tasdev_all_chroms_15samps_$cX.vcf 


	## --------------------------------
	## Apply final set of filters: data missingness, MAF
	## (not applying depth filters because they won't be applied to downsampled files)

	vcftools --vcf final_tasdev_all_chroms_15samps_$cX.vcf \
	--recode --recode-INFO-all \
	--out final_filt_tasdev_all_chroms_15samps__$cX \
	--maf 0.05 \
	--max-missing 0.8 \
	--minQ 20

done


























