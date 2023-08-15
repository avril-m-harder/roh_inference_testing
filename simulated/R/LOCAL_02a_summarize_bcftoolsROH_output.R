## Nomenclature note: 'true' ROHs are ROHs that were identified using sample FASTA files output from SLiM.
## 'called' ROHs are ROHs identified by bcftools/ROH from simulated reads aligned to the reference FASTA
## using various population sizes and coverage levels.

library(scales)

setwd('/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/bcftoolsroh_output/')

## get filenames of GT-only and PL results
fns <- list.files()
gt.fns <- fns[grep('_gt_', fns)]
pl.fns <- fns[grep('_pl_', fns)]

##### 1. Read in bcftools/ROH results, summarize, and write files with summary stats and ROH coordinates (can read in at #3) #####
## loop over GT-only results to summarize ROH call data
OUT <- NULL ## will store coordinates of all identified ROHs for all pop size / covg / sample
for(f in gt.fns){
  print(f)
  splt.f <- do.call(rbind, strsplit(f, split='_'))
  ## get pop size and coverage from filename
  pop.size <- splt.f[1,3]
  covg <- gsub('x', '', splt.f[1,5])
  ## read in data
  dat <- read.table(f)
  dat <- dat[,-1]
  colnames(dat) <- c('id','chrom','start','end','length','num.snps','quality')
  ## check for match here
  print(paste0('pop size = ',pop.size,' - number of samples = ',length(unique(dat$id))))
  for(i in unique(dat$id)){
    sub <- dat[dat$id == i,]
    ## save coordinates of all identified ROHs in every sample
    sub3 <- sub[,c(1,3:7)]
    sub3$pop.size <- pop.size
    sub3$covg <- covg
    sub3$id <- gsub('i','',sub3$id)
    sub3 <- sub3[,c(1,7,8,2:6)]
    OUT <- rbind(OUT, sub3)
  }
}
## write results
write.table(OUT, '/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/bcftoolsroh_results/bcftoolsROH_GT_all_coordinates.txt', 
            row.names=FALSE, quote = FALSE, sep = '\t')

## loop over PL results to summarize ROH call data
OUT <- NULL ## will store coordinates of all identified ROHs for all pop size / covg / sample
for(f in pl.fns){
  print(f)
  splt.f <- do.call(rbind, strsplit(f, split='_'))
  ## get pop size and coverage from filename
  pop.size <- splt.f[1,3]
  covg <- gsub('x', '', splt.f[1,5])
  ## read in data
  dat <- read.table(f)
  dat <- dat[,-1]
  colnames(dat) <- c('id','chrom','start','end','length','num.snps','quality')
  ## check for match here
  print(paste0('pop size = ',pop.size,' - number of samples = ',length(unique(dat$id))))
  for(i in unique(dat$id)){
    sub <- dat[dat$id == i,]
    ## save coordinates of all identified ROHs in every sample
    sub3 <- sub[,c(1,3:7)]
    sub3$pop.size <- pop.size
    sub3$covg <- covg
    sub3$id <- gsub('i','',sub3$id)
    sub3 <- sub3[,c(1,7,8,2:6)]
    OUT <- rbind(OUT, sub3)
  }
}
## write results
write.table(OUT, '/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/bcftoolsroh_results/bcftoolsROH_PL_all_coordinates.txt', 
          row.names=FALSE, quote = FALSE, sep = '\t')
