library(readr)
library(data.table)

correct_vcf <- function(vcf_file) {
      
meta <- readr::read_delim(vcf_file, delim = "\t", n_max = 5, 
                          col_names = FALSE, show_col_types = FALSE)[[1]]
names <- strsplit(vcf_file, split = '/')[[1]]
comps <- strsplit(names[length(names)], split = '_')
new.fn <- paste0('/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/',demo_mod,'/vcfs/final',comps[[1]][1],'_',demo_mod,'_',comps[[1]][2])

# read genotypes, data.table does this correctly
# i.e. jumping over the first 5 lines
gt <- data.table::fread(vcf_file, skip = 5)
# correct ref and alt columns
gt$REF <- 0
gt$ALT <- 1
# rename samples, can retrieve pedigree #s later if necessary
colnames(gt)[10:ncol(gt)] <- paste0('i',c(1:50))
# write out meta 
readr::write_lines(meta, file = new.fn)
# append genotypes to make file complete
data.table::fwrite(gt, file = new.fn, append = TRUE, 
       col.names = TRUE, sep = "\t")

}
