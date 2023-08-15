library(scales)

setwd('/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/plink_output_round2/')
fns <- list.files()

##### Read in PLINK results and summarize #####
## from shell script
# ## define parameter settings: round 1
# phwh=c(0, 1, 2)                        # Values for -homozyg-window-het
# phwm=c(2, 5, 50)                       # Values for -homozyg-window-missing
# phws=c(50, 100, 1000)                  # Values for -homozyg-window-snp
# phzd=c(50)                             # Values for -homozyg-density
# phzg=c(500, 1000)                      # Values for -homozyg-gap
# phwt=c(0.01, 0.05, 0.1)                # Values for -homozyg-window-threshold
# phzs=c(10, 100, 1000)                  # Values for -homozyg-snp
# phzk=c(100)                            # Values for -homozyg-kb
# pop.size <- c(100)                     # Values for sample size (EDITED)
# covg <- c('05','10','15','30','50')               # Values for mean coverage

## define parameter settings: round 2
phwh=c(1)                              # Values for -homozyg-window-het
phwm=c(5)                              # Values for -homozyg-window-missing
phws=c(40, 50, 60, 70, 80, 90, 100)    # Values for -homozyg-window-snp
phzd=c(50)                             # Values for -homozyg-density
phzg=c(1000)                           # Values for -homozyg-gap
phwt=c(0.05)                           # Values for -homozyg-window-threshold
phzs=c(10, 20, 40, 60, 80, 100)        # Values for -homozyg-snp
phzk=c(100)                            # Values for -homozyg-kb
pop.size <- c(100)                     # Values for sample size (EDITED)
covg <- c('05','10','15','30','50')    # Values for mean coverage

## Based on table (edited July 5, 2022), should be 2430 combos
nrow(expand.grid(pop.size, covg, phwh, phwm, phws, phzd, phzg, phwt, phzs, phzk))
#
##### Not using: Automated parameter setting retrieval #####
## For some reason, there were some extra leftover output files in the output directory. Manually setting the
## parameter settings of interest will weed these out.

# ## From list of file names, get lists of all parameter settings used in set of results files to be plotted
# splt.names <- do.call(rbind, strsplit(fns, split='_', fixed=TRUE))
# pop.ns <- do.call(rbind, strsplit(fns, split='pop_', fixed=TRUE))[,2]
# table(do.call(rbind, strsplit(pop.ns, split='_', fixed=TRUE))[,1])
# # pop.size <- unique(do.call(rbind, strsplit(pop.ns, split='_', fixed=TRUE))[,1])
# pop.size <- 100
# 
# covg.ns <- do.call(rbind, strsplit(fns, split='cvg_', fixed=TRUE))[,2]
# table(do.call(rbind, strsplit(covg.ns, split='_', fixed=TRUE))[,1])
# covg <- gsub('x', '', unique(do.call(rbind, strsplit(covg.ns, split='_', fixed=TRUE))[,1]))
# 
# phwh.ns <- do.call(rbind, strsplit(fns, split='phwh_', fixed=TRUE))[,2]
# table(do.call(rbind, strsplit(phwh.ns, split='_', fixed=TRUE))[,1])
# phwh <- unique(do.call(rbind, strsplit(phwh.ns, split='_', fixed=TRUE))[,1])
# 
# phwm.ns <- do.call(rbind, strsplit(fns, split='phwm_', fixed=TRUE))[,2]
# phwm <- unique(do.call(rbind, strsplit(phwm.ns, split='_', fixed=TRUE))[,1])
# 
# phws.ns <- do.call(rbind, strsplit(fns, split='phws_', fixed=TRUE))[,2]
# phws <- unique(do.call(rbind, strsplit(phws.ns, split='_', fixed=TRUE))[,1])
# 
# phzd.ns <- do.call(rbind, strsplit(fns, split='phzd_', fixed=TRUE))[,2]
# phzd <- unique(do.call(rbind, strsplit(phzd.ns, split='_', fixed=TRUE))[,1])
# 
# phzg.ns <- do.call(rbind, strsplit(fns, split='phzg_', fixed=TRUE))[,2]
# phzg <- unique(do.call(rbind, strsplit(phzg.ns, split='_', fixed=TRUE))[,1])
# 
# phwt.ns <- do.call(rbind, strsplit(fns, split='phwt_', fixed=TRUE))[,2]
# phwt <- unique(do.call(rbind, strsplit(phwt.ns, split='_', fixed=TRUE))[,1])
# 
# phzs.ns <- do.call(rbind, strsplit(fns, split='phzs_', fixed=TRUE))[,2]
# phzs <- unique(do.call(rbind, strsplit(phzs.ns, split='_', fixed=TRUE))[,1])
# 
# phzk.ns <- do.call(rbind, strsplit(fns, split='phzk_', fixed=TRUE))[,2]
# phzk <- unique(do.call(rbind, strsplit(phzk.ns, split='.', fixed=TRUE))[,1])




##### Summarize all coordinates in 1 file #####
z <- 1
OUT <- NULL
# OUT1 <- NULL
for(pop in pop.size){ 
  for(cov in covg){ 
    for(a in phwh){ 
      for(b in phwm){ 
        for(c in phws){ 
          for(d in phzd){
            for(e in phzg){
              for(f in phwt){
                for(g in phzs){
                  for(h in phzk){
                    fn <- paste0('sample_pop_',pop,'_cvg_',cov,'x_plink_roh_phwh_',a,'_phwm_',b,'_phws_',c,'_phzd_',d,'_phzg_',e,'_phwt_',f,'_phzs_',g,'_phzk_',h,'.hom')
                    print(z/7290)
                    z <- z+1
                    dat <- read.table(fn, header=TRUE)
                    if(nrow(dat) == 0){
                      # print(paste0('MISSING FILE: ',fn))
                      save <- c(pop, cov, a, b, c, d, e, f, g, h)
                      OUT <- rbind(OUT, save)
                      # OUT1 <- rbind(OUT1, fn)
                    } else{
                      dat$FID <- gsub('i', '', dat$FID)
                      dat <- dat[,c(1,7,8,10)]
                      dat$phwh <- a
                      dat$phwm <- b
                      dat$phws <- c
                      dat$phzd <- d
                      dat$phzg <- e
                      dat$phwt <- f
                      dat$phzs <- g
                      dat$phzk <- h
                      dat$pop.size <- pop
                      dat$covg <- cov
                      write.table(dat, '/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/plink_results_round2/PLINK_all_coordinates.txt',
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
}

# for(c in 1:ncol(OUT)){
#   print(table(OUT[,c]))
# }

#### After exploring the issue of missing results, 1280 runs of PLINK (corresponding to 1280 missing output files) found 0 ROH.
#### Confirmed by grep "found 0 ROH" on all of the slurm*out files and totaling up the number of matches.
#### 6010 combinations called >0 ROH.
#### See /Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/22_05_19_sorting_missing_plink_output_files/

## write table of combinations that yielded "found 0 ROH" in log files
colnames(OUT) <- c('pop.size', 'cov', 'phwh', 'phwm', 'phws', 'phzd', 'phzg', 'phwt', 'phzs', 'phzk')
write.csv(OUT, '../plink_results_round2/found_zero_rohs_combos.csv', row.names = FALSE)
