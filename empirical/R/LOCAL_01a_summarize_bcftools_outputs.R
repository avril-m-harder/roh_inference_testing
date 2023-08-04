`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/data/')

##### Read in ROH calling results and chrom data #####
## tas dev chrom data
chroms <- read.table('mSarHar1.11_autosomes.txt', sep = '\t', header = TRUE)
chroms <- chroms[,c(5,9)]
colnames(chroms) <- c('chrom','length')
tot.len <- sum(chroms$length)

## set categories
meths <- c('GT','PL')
covgs <- cbind(c(5, 10, 15, 30, 50), c('5X','10X','15X','30X','fullcovg'))
ests <- c('default','vtrained')

for(m in meths){
  OUT <- NULL
  for(c in 1:nrow(covgs)){
    for(e in ests){
      dat <- read.table(paste0('tasdev_',m,'only_',covgs[c,2],'_',e,'_RG_ONLY.txt'),
                        header = FALSE, sep = '\t')
      colnames(dat) <- c('RG','id','chrom','start','end','length','n.snps','quality')
      dat <- dat[dat$length >= 100e3,] ## applying 100-kb ROH length filter
      dat$covg <- covgs[c,1]
      dat$ests <- e
      OUT <- rbind(OUT, dat)
    }
  }
  write.csv(OUT[,c(2:ncol(OUT))], paste0('tasdev_',m,'_allcovg_levels.csv'), row.names = FALSE)
}
