setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/data/')

##### 1. Read in data and summary files #####
## Read in sample and chrom keys
samp.key <- read.csv('../../../empirical_qc/sample_id_key.csv')
chrom.key <- read.csv('../../../empirical_qc/chrom_key.csv')

## tas dev chrom data
chroms <- read.table('mSarHar1.11_autosomes.txt', sep = '\t', header = TRUE)
chroms <- chroms[,c(5,9)]
colnames(chroms) <- c('chrom','length')
tot.len <- sum(chroms$length)

## BCFtools/ROH results
gt.res <- read.csv('tasdev_GT_allcovg_levels.csv')
gt.res <- gt.res[gt.res$ests == 'vtrained',]
gt.res <- merge(gt.res, samp.key, by.x = 'id', by.y = 'sra.id')
colnames(gt.res)[c(1,ncol(gt.res))] <- c('sra.id','id')
pl.res <- read.csv('tasdev_PL_allcovg_levels.csv')
pl.res <- pl.res[pl.res$ests == 'vtrained',]
pl.res <- merge(pl.res, samp.key, by.x = 'id', by.y = 'sra.id')
colnames(pl.res)[c(1,ncol(pl.res))] <- c('sra.id','id')

## PLINK results
plink.res <- read.csv('PLINK_all_coordinates_DEFAULT_SETTINGS.csv')
plink.res$length <- plink.res$end - plink.res$start + 1

###
OUT <- NULL
for(i in unique(gt.res$id)){
  for(c in unique(gt.res$covg)){
    sub.gt <- gt.res[gt.res$id == i & gt.res$covg == c,]
    sub.pl <- pl.res[pl.res$id == i & pl.res$covg == c,]
    sub.plink <- plink.res[plink.res$num.id == i & plink.res$covg == c,]
    
    save <- c(i, c, sum(sub.gt$length)/tot.len, sum(sub.pl$length)/tot.len, sum(sub.plink$length)/tot.len)
    OUT <- rbind(OUT, save)
  }
}
colnames(OUT) <- c('num.id','covg','gt.froh','pl.froh','plink.froh')
write.csv(OUT, 'individual_froh_stats.csv', row.names = FALSE)
