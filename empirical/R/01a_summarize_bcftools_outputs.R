`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Desktop/roh_param_project/empirical_plotting_summarizing/')

##### Read in ROH calling results and chrom data #####
## tas dev chrom data
chroms <- read.table('tasdev_assembly_data/mSarHar1.11_autosomes.txt', sep = '\t', header = TRUE)
chroms <- chroms[,c(5,9)]
colnames(chroms) <- c('chrom','length')
tot.len <- sum(chroms$length)

## BCFtools/ROH
## GT
## Full coverage
bcf.gt.full <- read.table('bcftools_data/tasdev_GTonly_fullcovg_RG_ONLY.txt',
                         header = FALSE, sep = '\t')
colnames(bcf.gt.full) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.gt.full <- bcf.gt.full[bcf.gt.full$length >= 100000,] ## applying 100kb filter, as would normally be done
bcf.gt.full$covg <- 50
## 5X
bcf.gt.5X <- read.table('bcftools_data/tasdev_GTonly_5X_RG_ONLY.txt',
                          header = FALSE, sep = '\t')
colnames(bcf.gt.5X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.gt.5X <- bcf.gt.5X[bcf.gt.5X$length >= 100000,]
bcf.gt.5X$covg <- 5
## 10X
bcf.gt.10X <- read.table('bcftools_data/tasdev_GTonly_10X_RG_ONLY.txt',
                        header = FALSE, sep = '\t')
colnames(bcf.gt.10X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.gt.10X <- bcf.gt.10X[bcf.gt.10X$length >= 100000,]
bcf.gt.10X$covg <- 10
## 15X
bcf.gt.15X <- read.table('bcftools_data/tasdev_GTonly_15X_RG_ONLY.txt',
                        header = FALSE, sep = '\t')
colnames(bcf.gt.15X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.gt.15X <- bcf.gt.15X[bcf.gt.15X$length >= 100000,]
bcf.gt.15X$covg <- 15
## 30X
bcf.gt.30X <- read.table('bcftools_data/tasdev_GTonly_30X_RG_ONLY.txt',
                        header = FALSE, sep = '\t')
colnames(bcf.gt.30X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.gt.30X <- bcf.gt.30X[bcf.gt.30X$length >= 100000,]
bcf.gt.30X$covg <- 30
## combine all coverage levels
bcf.gt.res <- rbind(bcf.gt.full, bcf.gt.5X, bcf.gt.10X, bcf.gt.15X, bcf.gt.30X)
write.csv(bcf.gt.res[,c(2:ncol(bcf.gt.res))], 'bcftools_data/tasdev_GT_allcovg_levels.csv', row.names = FALSE)

## PL
## Full coverage
bcf.pl.full <- read.table('bcftools_data/tasdev_PLonly_fullcovgL_RG_ONLY.txt',
                         header=FALSE, sep = '\t')
colnames(bcf.pl.full) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.pl.full <- bcf.pl.full[bcf.pl.full$length >= 100000,] ## applying 100kb filter, as would normally be done
bcf.pl.full$covg <- 50
## 5X
bcf.pl.5X <- read.table('bcftools_data/tasdev_PLonly_5X_RG_ONLY.txt',
                        header = FALSE, sep = '\t')
colnames(bcf.pl.5X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.pl.5X <- bcf.pl.5X[bcf.pl.5X$length >= 100000,]
bcf.pl.5X$covg <- 5
## 10X
bcf.pl.10X <- read.table('bcftools_data/tasdev_PLonly_10X_RG_ONLY.txt',
                        header = FALSE, sep = '\t')
colnames(bcf.pl.10X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.pl.10X <- bcf.pl.10X[bcf.pl.10X$length >= 100000,]
bcf.pl.10X$covg <- 10
## 15X
bcf.pl.15X <- read.table('bcftools_data/tasdev_PLonly_15X_RG_ONLY.txt',
                        header = FALSE, sep = '\t')
colnames(bcf.pl.15X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.pl.15X <- bcf.pl.15X[bcf.pl.15X$length >= 100000,]
bcf.pl.15X$covg <- 15
## 30X
bcf.pl.30X <- read.table('bcftools_data/tasdev_PLonly_30X_RG_ONLY.txt',
                        header = FALSE, sep = '\t')
colnames(bcf.pl.30X) <- c('RG','id','chrom','start','end','length','n.snps','quality')
bcf.pl.30X <- bcf.pl.30X[bcf.pl.30X$length >= 100000,]
bcf.pl.30X$covg <- 30
## combine all coverage levels
bcf.pl.res <- rbind(bcf.pl.full, bcf.pl.5X, bcf.pl.10X, bcf.pl.15X, bcf.pl.30X)
write.csv(bcf.pl.res[,c(2:ncol(bcf.pl.res))], 'bcftools_data/tasdev_PL_allcovg_levels.csv', row.names = FALSE)