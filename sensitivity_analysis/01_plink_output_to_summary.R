library(scales)

## set round and data set
set <- 'empirical' ### empirical    simulated
# rds <- c('round1','round2a','round2b','round2c','round2d','round2e')
rds <- c('round2a','round2b','round2c','round2d','round2e')

for(rd in rds){
  setwd(paste0('/Users/Avril/Desktop/roh_param_project/',set,'_plotting_summarizing/plink_data_',rd,'/'))
  fns <- list.files()
  
  ##### Summarize all coordinates in 1 file #####
  z <- 1
  for(fn in fns){
    print(z/length(fns))
    ## read in data
    dat <- read.table(fn, header=TRUE)
    
    if(nrow(dat) > 0){
      ## split file name to get parameter settings
      covg.ns <- do.call(rbind, strsplit(fn, split='cvg_', fixed=TRUE))[,2]
      cov <- gsub('x', '', unique(do.call(rbind, strsplit(covg.ns, split='_', fixed=TRUE))[,1]))
    
      phwh.ns <- do.call(rbind, strsplit(fn, split='phwh_', fixed=TRUE))[,2]
      a <- unique(do.call(rbind, strsplit(phwh.ns, split='_', fixed=TRUE))[,1])
    
      phwm.ns <- do.call(rbind, strsplit(fn, split='phwm_', fixed=TRUE))[,2]
      b <- unique(do.call(rbind, strsplit(phwm.ns, split='_', fixed=TRUE))[,1])
    
      phws.ns <- do.call(rbind, strsplit(fn, split='phws_', fixed=TRUE))[,2]
      c <- unique(do.call(rbind, strsplit(phws.ns, split='_', fixed=TRUE))[,1])
    
      phzd.ns <- do.call(rbind, strsplit(fn, split='phzd_', fixed=TRUE))[,2]
      d <- unique(do.call(rbind, strsplit(phzd.ns, split='_', fixed=TRUE))[,1])
    
      phzg.ns <- do.call(rbind, strsplit(fn, split='phzg_', fixed=TRUE))[,2]
      e <- unique(do.call(rbind, strsplit(phzg.ns, split='_', fixed=TRUE))[,1])
    
      phwt.ns <- do.call(rbind, strsplit(fn, split='phwt_', fixed=TRUE))[,2]
      f <- unique(do.call(rbind, strsplit(phwt.ns, split='_', fixed=TRUE))[,1])
    
      phzs.ns <- do.call(rbind, strsplit(fn, split='phzs_', fixed=TRUE))[,2]
      g <- unique(do.call(rbind, strsplit(phzs.ns, split='_', fixed=TRUE))[,1])
    
      phzk.ns <- do.call(rbind, strsplit(fn, split='phzk_', fixed=TRUE))[,2]
      h <- unique(do.call(rbind, strsplit(phzk.ns, split='_', fixed=TRUE))[,1])
      h <- gsub('.hom', '' ,h)
      
      ## prep data to keep and write
      dat$FID <- gsub('i', '', dat$FID)
      if(set == 'simulated'){
        dat <- dat[,c(1,7,8,10)]
      }
      if(set == 'empirical'){
        dat <- dat[,c(1,4,7,8,10)]
      }
      dat$phwh <- a
      dat$phwm <- b
      dat$phws <- c
      dat$phzd <- d
      dat$phzg <- e
      dat$phwt <- f
      dat$phzs <- g
      dat$phzk <- h
      dat$covg <- cov
      
      write.table(dat, paste0('/Users/Avril/Desktop/roh_param_project/',set,'_plotting_summarizing/plink_results_',rd,'/PLINK_all_coordinates_',rd,'.txt'),
                  append = TRUE, quote = FALSE, row.names = FALSE, sep='\t', col.names = FALSE)
    } 
    z <- z+1
  }
  
  
  ##### Summarize individual f(ROH) data #####
  if(set == 'simulated'){
    true.rohs <- read.table('../true_indiv_data/true_roh_coords.txt')
    colnames(true.rohs) <- c('id','start','end','length')
    true.rohs <- true.rohs[true.rohs$length >= 100000,] ### only keeping true ROHs >= 100 kb because that's all we're evaluating the ability to call
    ## create unique ROH ID for linking true ROHs to called ROHs
    ns <- c(1:nrow(true.rohs))
    true.rohs$true.roh.id <- ns
    chrom.len <- 30e6
    
    plink.out <- read.table(paste0('../plink_results_',rd,'/PLINK_all_coordinates_',rd,'.txt'))
    colnames(plink.out) <- c('id','start','end','n.snps','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','covg')
    
    OUT <- NULL
    write.table(OUT, paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.txt'),
                sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
    for(cov in unique(plink.out$covg)){
      for(a in unique(plink.out$phwh)){
        for(b in unique(plink.out$phwm)){
          for(c in unique(plink.out$phws)){
            for(d in unique(plink.out$phzd)){
              for(e in unique(plink.out$phzg)){
                for(f in unique(plink.out$phwt)){
                  for(g in unique(plink.out$phzs)){
                    for(h in unique(plink.out$phzk)){
                      sub <- plink.out[plink.out$covg == cov & plink.out$phwh == a &
                                       plink.out$phwm == b & plink.out$phws == c & plink.out$phzd == d &
                                       plink.out$phzg == e & plink.out$phwt == f & plink.out$phzs == g &
                                       plink.out$phzk == h,]
                      if(nrow(sub) > 0){
                        for(i in unique(sub$id)){
                          temp <- sub[sub$id == i,]
                          temp$length <- temp$end - temp$start + 1
                          true.froh <- sum(true.rohs[true.rohs$id == i, 'length'])/30e6
                          call.froh <- sum(temp$length)/30e6
                          save <- c(i, cov, a, b, c, d, e, f, g, h, true.froh, call.froh)
                          OUT <- rbind(OUT, save)
                        }
                      }
                      write.table(OUT, paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.txt'),
                                  sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
                      OUT <- NULL
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  
    OUT <- read.table(paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.txt'))
    colnames(OUT) <- c('id','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','true.froh','call.froh')
    froh.stats <- as.data.frame(OUT)
    froh.stats$abs.diff <- abs(froh.stats$call.froh - froh.stats$true.froh)
    froh.stats$group <- paste0(froh.stats$pop.size,'-',froh.stats$covg,'-',froh.stats$phwh,'-',froh.stats$phwm,'-',froh.stats$phws,'-',froh.stats$phzd,'-',froh.stats$phzg,'-',froh.stats$phwt,'-',froh.stats$phzs,'-',froh.stats$phzk)
    froh.stats <- froh.stats[,c('id','true.froh','call.froh','abs.diff','group','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk')]
    write.csv(froh.stats, paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.csv'), row.names = FALSE)
    write.csv(froh.stats, paste0('/Users/Avril/Desktop/roh_param_project/iterative_sa_steps/',set,'/',rd,'/individual_froh_results_',rd,'.csv'), row.names = FALSE)
    froh.stats <- read.csv(paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.csv'))
  }
  if(set == 'empirical'){
    ## tas dev chrom data
    chroms <- read.table('../tasdev_assembly_data/mSarHar1.11_autosomes.txt', sep = '\t', header = TRUE)
    chroms <- chroms[,c(5,9)]
    colnames(chroms) <- c('chrom','length')
    tot.len <- sum(chroms$length)
    
    plink.out <- read.table(paste0('../plink_results_',rd,'/PLINK_all_coordinates_',rd,'.txt'))
    colnames(plink.out) <- c('id','chrom','start','end','n.snps','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','covg')
    covs <- cbind(c('5X','10X','15X','30X','fullcovg'), c(5, 10, 15, 30, 50))
    colnames(covs) <- c('char','num')
    plink.out <- merge(plink.out, covs, by.x = 'covg', by.y = 'char')
    plink.out <- plink.out[,-1]
    colnames(plink.out)[ncol(plink.out)] <- 'covg'
    
    OUT <- NULL
    write.table(OUT, paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.txt'),
                sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
    for(cov in unique(plink.out$covg)){
      for(a in unique(plink.out$phwh)){
        for(b in unique(plink.out$phwm)){
          for(c in unique(plink.out$phws)){
            for(d in unique(plink.out$phzd)){
              for(e in unique(plink.out$phzg)){
                for(f in unique(plink.out$phwt)){
                  for(g in unique(plink.out$phzs)){
                    for(h in unique(plink.out$phzk)){
                      sub <- plink.out[plink.out$covg == cov & plink.out$phwh == a &
                                         plink.out$phwm == b & plink.out$phws == c & plink.out$phzd == d &
                                         plink.out$phzg == e & plink.out$phwt == f & plink.out$phzs == g &
                                         plink.out$phzk == h,]
                      if(nrow(sub) > 0){
                        for(i in unique(sub$id)){
                          temp <- sub[sub$id == i,]
                          temp$length <- temp$end - temp$start + 1
                          call.froh <- sum(temp$length)/tot.len
                          save <- c(i, cov, a, b, c, d, e, f, g, h, call.froh)
                          OUT <- rbind(OUT, save)
                        }
                      }
                      write.table(OUT, paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.txt'),
                                  sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
                      OUT <- NULL
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    OUT <- read.table(paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.txt'))
    colnames(OUT) <- c('id','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','call.froh')
    froh.stats <- as.data.frame(OUT)
    froh.stats <- froh.stats[,c('id','call.froh','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk')]
    write.csv(froh.stats, paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.csv'), row.names = FALSE)
    write.csv(froh.stats, paste0('/Users/Avril/Desktop/roh_param_project/iterative_sa_steps/',set,'/',rd,'/individual_froh_results_',rd,'.csv'), row.names = FALSE)
    froh.stats <- read.csv(paste0('../plink_results_',rd,'/individual_froh_results_',rd,'.csv'))
  }
}