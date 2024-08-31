# ROH Inference Testing

<p style="font-size: 20px;">
Scripts and Files Associated with the Manuscript: "Detectability of Runs of Homozygosity is Influenced by Analysis Parameters and Population-Specific Demographic History"
</p>

This repository contains scripts and data files related to the analysis conducted in the manuscript currently available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.09.29.510155v1) and submitted for peer-review.

---

## Empirical Data - Bash Scripts

- **01a_download_and_qc.sh**: Downloads read files from SRA, performs quality checks with FastQC, and cleans data with TrimGalore.
- **01b_assem_indexing.sh**: Downloads and indexes the reference assembly.
- **02a_read_mapping_indivsamples.sh**: Maps reads to the reference genome using BWA-MEM.
- **02b_downsampling.sh**: Downsamples full coverage sets to 5X, 10X, 15X, and 30X using SAMtools.
- **03a_preBQSR_snp_calling_fullcovg.sh**: Prepares BAM files for HaplotypeCaller and runs the process for full coverage data.
- **03a_preBQSR_snp_calling.sh**: Similar to the above but for downsampled coverage levels (5X-30X).
- **03b_noBQSR_genomicsdb_chrom_batches_fullcovg.sh**: Genotypes GVCFs and applies hard filtering in GATK and VCFtools for full coverage data.
- **03b_noBQSR_genomicsdb_chrom_batches.sh**: Same as the above but for downsampled coverage levels (5X-30X).
- **04a_bcftoolsROH.sh**: Runs BCFtools/RoH in both Likelihoods and Genotypes modes.
- **05b_plink_empirical_round1.sh**: Runs PLINK with a defined set of parameter settings. Parameters were adjusted in subsequent iterations.
- **depth_check.sh**: Checks coverage in BAM files using SAMtools.
- **indiv_samp_script_gen.sh**: Generates and submits individual sample scripts for SNP calling.
- **qualimap.sh**: Runs Qualimap at various BAM processing stages.
- **tasdev_vcf_sample_depths.sh**: Checks read depths at final sets of filtered SNPs.

## Empirical Data - R Scripts

- **01a_summarize_bcftools_outputs.R**: Summarizes BCFtools outputs (ROH coordinates) across all coverage levels for input into downstream analyses.
- **01b_compare_viterbi_default.R**: Compares BCFtools ROH calls with and without the Viterbi algorithm.
- **01c_summarize_plink_output.R**: Summarizes PLINK results (ROH coordinates) across all coverage levels for input into downstream analyses.
- **02_indiv_fROH_calcs.R**: Calculates individual *F*<sub>ROH</sub> values from all methods and coverage levels for downstream analysis.
- **03_manuscript_figures_and_statistics_empirical.R**: Conducts all downstream analyses and generates corresponding figures for the manuscript using results from all methods and coverage levels.
- **04_visualizing_het.R**: Calculates and visualizes heterozygosity in windows alongside true and called ROHs.

## Simulated Data

### Bash Scripts

**Note:** SLiM and bash scripts for simulated data generation and analyses are available in [kennethb22/roh_parameter_project_kk](https://github.com/kennethb22/roh_parameter_project_kk). Original scripts were written by Kenneth Kirksey; modified versions are included in this repository as listed bellow.

- **01_slim**
    - **01_process_slim_output.sh.sh**: Processes the simulated data from SLiM.
- **02_read_sim_and_qc**
    - **02a_run_art_w_make.sh**: Runs ART to simulate reads.
    - **02b_concat_fastq_w_make.sh**: Concatenates FASTQ files.
    - **02c_fastqc_w_make.sh**: Evaluates reads quality with FastQC.
- **03_read_align**
    - **03_run_bwa_w_make.sh**: Aligns reads to the reference genome.
- **04_downsample**
    - **04_run_downsample_w_make.sh**: Downsamples BAM files to 5X, 10X, 15X, and 30X.
- **05_var_call**
    - **05a_add_read_groups.sh**: Adds read groups to BAM files.
    - **05b_haplotype_calling.sh**: Calls SNPs with HaplotypeCaller.
    - **05c_combine_and_genotype_gvcfs_w_make.sh**: Combines GVCFs and genotypes.
    - **05d_filter_variants_w_make.sh**: Filters variants.
- **06_roh_id**
    - **06a_run_bcftoolsROH_w_make.sh**: Runs BCFtools ROH.
    - **06a_run_bcftoolsROH_w_sample_exclusion.sh**: Runs BCFtools ROH with sample exclusion.
    - **06b_run_plink_w_make.sh**: Runs PLINK.

### R Scripts

- **01_getting_true_roh_coords.R**: Summarizes the start and end coordinates of true ROHs.
- **02a_summarize_bcftoolsROH_output.R**: Summarizes BCFtools outputs (ROH coordinates) across all coverage levels for input into downstream analyses.
- **02b_summarize_plink_output.R**: Summarizes PLINK results (ROH coordinates) across all coverage levels for input into downstream analyses.
- **03a_identifying_overlap_bcftoolsROH.R**: Identifies overlapping regions between true and called ROH coordinates.
- **03b_roh_param_analyses_preSA_plink.R**: Summarizes PLINK results and identifies overlap between true and called ROHs.
- **04a_bcftools_evaluating_viterbi_effects.R**: Evaluates the effects of the Viterbi algorithm on BCFtools ROH calls.
- **05_manuscript_figures_and_statistics_simulated.R**: Conducts all downstream analyses and generates corresponding figures for the manuscript using results from all methods and coverage levels.
- **06_visualizing_het.R**: Calculates and visualizes heterozygosity in windows alongside true and called ROHs


## Sensitivity Analysis - R Scripts

- **01_plink_output_to_summary.R**: Reads in results from one iteration of PLINK setting combinations and writes input for sensitivity analysis.
- **02_sensitivity_analysis.R**: Conducts sensitivity analysis and produces diagnostic plots.

## Files

- **15_samples.txt**: Accession numbers for SRA read downloads.
- **individual_froh_results_SA.csv**: Example input file for `02_sensitivity_analysis.R`.

---