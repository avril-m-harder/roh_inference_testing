# roh_inference_testing

### Scripts and files associated with ROH inference testing manuscript, currently available on [bioRxiv](link goes in these parentheses) and submitted for peer-review.

#### Empirical data - bash scripts
* *01a_download_and_qc.sh*: Download read files from SRA, check quality with FastQC, and clean with TrimGalore
* *01b_assem_indexing.sh*: Download and index reference assembly
* *02a_read_mapping_indivsamples.sh*: Map.... reads. With BWA-MEM
* *02b_downsampling.sh*: Downsample full coverage set to 5X, 10X, 15X, and 30X coverage levels with SAMtools
* *03a_preBQSR_snp_calling_fullcovg.sh*: Prep BAM files for HaplotypeCaller and run 
* *03a_preBQSR_snp_calling.sh*: Same as 03a_preBQSR_snp_calling_fullcovg.sh but for 5X-30X coverage levels
* *03b_noBQSR_genomicsdb_chrom_batches_fullcovg.sh*: Genotype GVCFs and apply hard filtering in GATK and VCFtools
* *03b_noBQSR_genomicsdb_chrom_batches.sh*: Same as 03b_noBQSR_genomicsdb_chrom_batches_fullcovg.sh but for 5X-30X coverage levels
* *04a_bcftoolsROH.sh*: Run BCFtools/RoH in Likelihoods and Genotypes modes
* *05b_plink_empirical_round1.sh*: Run PLINK with a defined set of parameter settings. PLINK parameters section was edited in subsequent iterations
* *depth_check.sh*: Check coverage in BAM files using samtools
* *indiv_samp_script_gen.sh*: Generate and submit 03a_preBQSR_snp_calling_\*.sh scripts for separate samples
* *qualimap.sh*: Run Qualimap at various BAM processing stages
* *tasdev_vcf_sample_depths.sh*: Check read depths at final sets of filtered SNPs
 
#### Empirical data - R scripts
* *01a_summarize_bcftools_outputs.R*: Reads in all of the output from BCFtools Genotypes and Likelihoods (i.e., ROH coordinates) across all coverage levels and summarizes into a single file for reading in with 03_manuscript_figures_and_statistics_empirical.R
* *01b_summarize_plink_output.R*: Reads in final set of PLINK results (i.e., ROH coordinates) for all coverage levels and summarizes for reading in with 03_manuscript_figures_and_statistics_empirical.R
* *02_indiv_fROH_calcs.R*: Reads in result summaries from all 3 methods and coverage levels and calculates individual *F*<sub>ROH</sub> values for reading in with 03_manuscript_figures_and_statistics_empirical.R
* *03_manuscript_figures_and_statistics_empirical.R*: Reads in results for all 3 methods and coverage levels and conducts all downstream analyses and produces all figures for the manuscript (plus a ton of.. bonus figures)

#### Simulated data scripts
* 

* 

* 

* 

* 

 


#### Files
* 15_samples.txt: Accession numbers for SRA read downloads.

* 
