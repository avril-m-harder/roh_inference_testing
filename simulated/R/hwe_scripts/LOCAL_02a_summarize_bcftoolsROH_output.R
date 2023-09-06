## Nomenclature note: 'true' ROHs are ROHs that were identified using sample FASTA files output from SLiM.
## 'called' ROHs are ROHs identified by bcftools/ROH from simulated reads aligned to the reference FASTA
## using various population sizes and coverage levels.

library(scales)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/bcftools_output/hwe_filtered/')

## get filenames of GT-only and PL results
fns <- list.files(pattern = 'RG_ONLY')

##### Read in bcftools/ROH results, summarize, and write files with summary stats and ROH coordinates #####
## loop over GT-only results to summarize ROH call data
OUT <- NULL ## will store coordinates of all identified ROHs for all method / demo scenario / covg / sample

called.roh.id.counter <- 1

for(f in fns){
  print(f)
  comps <- unlist(strsplit(f, split = '_'))
  demo <- comps[1]
  method <- comps[2]
  covg <- comps[3]
  hmm <- comps[6]
  ## read in data
  dat <- read.table(f)
  dat <- dat[,-1]
  colnames(dat) <- c('id','chrom','start','end','length','num.snps','quality')
  s.id <- called.roh.id.counter
  e.id <- called.roh.id.counter + nrow(dat) - 1
  dat$called.roh.id <- c(s.id:e.id)
  called.roh.id.counter <- called.roh.id.counter + nrow(dat)
  ## summarize and save data
  for(i in unique(dat$id)){
    sub <- dat[dat$id == i,]
    ## save coordinates of all identified ROHs in every sample
    sub3 <- sub[,c(1,3:8)]
    sub3$demo <- demo
    sub3$method <- method
    sub3$covg <- covg
    sub3$hmm <- hmm
    sub3 <- sub3[,c(1,8:11,2:7)]
    OUT <- rbind(OUT, sub3)
  }
}
write.table(OUT, 'bcftoolsROH_all_coordinates.txt', 
            row.names=FALSE, quote = FALSE, sep = '\t')
