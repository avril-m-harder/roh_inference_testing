library(scales)

rd <- 'round2b'

setwd(paste0('/Users/Avril/Desktop/roh_param_project/empirical_plotting_summarizing/plink_data_',rd,'/'))
fns <- list.files()

##### Read in sample and chrom keys #####
samp.key <- read.csv('../../empirical_qc/sample_id_key.csv')
chrom.key <- read.csv('../../empirical_qc/chrom_key.csv')

##### Read in PLINK results and summarize #####
## from shell script
## define parameter settings: round 1
# phwh=c(0, 1, 2)                               # Values for -homozyg-window-het
# phwm=c(2, 5, 50)                              # Values for -homozyg-window-missing
# phws=c(50, 100, 1000)                         # Values for -homozyg-window-snp
# phzd=c(50)                                    # Values for -homozyg-density
# phzg=c(500, 1000)                             # Values for -homozyg-gap
# phwt=c(0.01, 0.05, 0.1)                       # Values for -homozyg-window-threshold
# phzs=c(10, 100, 1000)                         # Values for -homozyg-snp
# phzk=c(100)                                   # Values for -homozyg-kb
# covg <- cbind(c('5X','10X','15X','30X','fullcovg'),
#               c(5,10,15,30,50))                # Values for mean coverage

## define parameter settings: round 2a
# phwh=c(1)                              # Values for -homozyg-window-het
# phwm=c(5)                              # Values for -homozyg-window-missing
# phws=c(50, 75, 100, 125)    # Values for -homozyg-window-snp
# phzd=c(50)                             # Values for -homozyg-density
# phzg=c(1000)                           # Values for -homozyg-gap
# phwt=c(0.05)                           # Values for -homozyg-window-threshold
# phzs=c(10, 20, 40, 60, 80, 100)        # Values for -homozyg-snp
# phzk=c(100)                            # Values for -homozyg-kb
# covg <- cbind(c('5X','10X','15X','30X','fullcovg'),
#               c(5,10,15,30,50))        # Values for mean coverage

## define parameter settings: round 2b
phwh=c(1)                              # Values for -homozyg-window-het
phwm=c(5)                              # Values for -homozyg-window-missing
phws=c(925, 950, 975, 1000)    # Values for -homozyg-window-snp
phzd=c(50)                             # Values for -homozyg-density
phzg=c(1000)                           # Values for -homozyg-gap
phwt=c(0.05)                           # Values for -homozyg-window-threshold
phzs=c(900, 920, 940, 960, 980, 1000)        # Values for -homozyg-snp
phzk=c(100)                            # Values for -homozyg-kb
covg <- cbind(c('5X','10X','15X','30X','fullcovg'),
              c(5,10,15,30,50))        # Values for mean coverage

## Based on table (edited July 5, 2022), should be 2430 combos for round 1
nrow(expand.grid(covg[,1], phwh, phwm, phws, phzd, phzg, phwt, phzs, phzk))


##### Summarize all coordinates in 1 file #####
z <- 1
OUT <- NULL
# OUT1 <- NULL
for(r in 1:nrow(covg)){ 
  for(a in phwh){ 
    for(b in phwm){ 
      for(c in phws){ 
        for(d in phzd){
          for(e in phzg){
            for(f in phwt){
              for(g in phzs){
                for(h in phzk){
                  fn <- paste0('tasdev_cvg_',covg[r,1],'_plink_roh_phwh_',a,'_phwm_',b,'_phws_',c,'_phzd_',d,'_phzg_',e,'_phwt_',f,'_phzs_',g,'_phzk_',h,'_',rd,'.hom')
                  print(z/nrow(expand.grid(covg[,1], phwh, phwm, phws, phzd, phzg, phwt, phzs, phzk)))
                  z <- z+1
                  dat <- read.table(fn, header=TRUE)
                  dat <- merge(dat, samp.key, by.x = 'FID', by.y = 'sra.id')
                  dat <- merge(dat, chrom.key, by.x = 'CHR', by.y = 'chrom.name')
                  if(nrow(dat) == 0){
                    # print(paste0('MISSING FILE: ',fn))
                    save <- c(covg[r,2], a, b, c, d, e, f, g, h)
                    OUT <- rbind(OUT, save)
                    # OUT1 <- rbind(OUT1, fn)
                  } else{
                    dat <- dat[,c(14,15,7,8,10)]
                    dat$phwh <- a
                    dat$phwm <- b
                    dat$phws <- c
                    dat$phzd <- d
                    dat$phzg <- e
                    dat$phwt <- f
                    dat$phzs <- g
                    dat$phzk <- h
                    dat$covg <- covg[r,2]
                    write.table(dat, paste0('/Users/Avril/Desktop/roh_param_project/empirical_plotting_summarizing/plink_results_',rd,'/PLINK_all_coordinates.txt'),
                                append = TRUE, quote = FALSE, row.names = FALSE, sep='\t', col.names = FALSE)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


# for(c in 1:ncol(OUT)){
#   print(table(OUT[,c]))
# }

## write table of combinations that yielded "found 0 ROH" in log files
# colnames(OUT) <- c('cov', 'phwh', 'phwm', 'phws', 'phzd', 'phzg', 'phwt', 'phzs', 'phzk')
# write.csv(OUT, '../plink_results_round1/found_zero_rohs_combos.csv', row.names = FALSE)
