library(scales)
library(grDevices)
`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/plink_output/')

## 
rnd <- 'hwe_round1'

##### Loop over demographic scenarios #####
true.fns <- list.files('../slim_true_data/', pattern = 'true_roh_coords')

for(true.fn in true.fns){
  demo <- unlist(strsplit(true.fn, split = '_'))[1]
  
  ## set up output files
  PLINK.OUT <- matrix(c('id','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk',
                        'len','true.len','true.roh.id','called.roh.id','n.snps','overlap.start','overlap.end'), nrow = 1)
  suppressWarnings(write.table(PLINK.OUT, paste0(rnd,'/',demo,'_PLINK_overlap_results.txt'),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE))
  
  OUT <- matrix(c('id','covg','phwh','phwm','phws','phzd','phzg',
                  'phwt','phzs','phzk','true.froh','call.froh'), nrow = 1)
  suppressWarnings(write.table(OUT, paste0(rnd,'/',demo,'_PLINK_individual_froh_results.txt'),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE))
  
  froh.stats <- matrix(c('id','true.froh','call.froh','abs.diff','combo','covg',
                         'phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk'), nrow = 1)
  suppressWarnings(write.csv(froh.stats, paste0(rnd,'/',demo,'_PLINK_individual_froh_results.csv'),
                             row.names = FALSE, col.names = FALSE))
  
  OUT1 <- matrix(c('mean.true.froh','mean.called.froh','sd.diff','group','covg',
                   'phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk'), nrow = 1)
  suppressWarnings(write.csv(OUT1, paste0(rnd,'/',demo,'_population_froh_results_SA.csv'), 
            row.names = FALSE, col.names = FALSE))
  
  froh.stats <- matrix(c('id','true.froh','call.froh','covg','phwh','phwm','phws',
                         'phzd','phzg','phwt','phzs','phzk'), nrow = 1)
  suppressWarnings(write.csv(froh.stats, paste0(rnd,'/',demo,'_individual_froh_results_SA.csv'), 
            row.names = FALSE, col.names = FALSE))
    
  ##### Read in known/true heterozygosity + ROH information #####
  true.rohs <- read.table(paste0('../slim_true_data/',true.fn))
  colnames(true.rohs) <- c('id','start','end','length','true.roh.id')
  true.rohs <- true.rohs[true.rohs$length >= 100000,] ### only keeping true ROHs >= 100 kb because that's all we're evaluating the ability to call
  chrom.len <- 30e6
  if(demo == 'decline'){
    true.rohs <- true.rohs[true.rohs$id %notin% c(33,48),]
  }
  
  ##### Read in ROH calling results #####
  ## PLINK (written from 02b_summarize_plink_output.R)
  plink.res <- read.table(paste0(rnd,'/',demo,'_PLINK_all_coordinates.txt'), header = TRUE)
  plink.res$id <- gsub('i','',plink.res$id)
  
  ## create combo variable to loop over
  plink.res$combo <- paste0(plink.res[,5],'-',
                            plink.res[,6],'-',
                            plink.res[,7],'-',
                            plink.res[,8],'-',
                            plink.res[,9],'-',
                            plink.res[,10],'-',
                            plink.res[,11],'-',
                            plink.res[,12],'-',
                            plink.res[,13])
  
  if(demo == 'decline'){
    plink.res <- plink.res[plink.res$id %notin% c(33,48),]
  }
  
  ##### Loop over param setting combos, writing ROHverlap results as it goes #####
  for(c in unique(plink.res$combo)){
    
    print(paste0(demo,' - ',c))
    
    ##### 2. ROH identification -- PLINK results #####
    ## Output: for each true ROH, identify all overlapping called ROHs. Each instance
    ## of overlap will occupy one line in the output matrix. The output matrix will contain:
    ## roh.id for true ROH, id/pop.size/covg/start/end for called ROH, and length of overlap.
    
    plink.combo <- plink.res[plink.res$combo == c,]
    
    ## For each true ROH, calculate...
    ## the proportion overlapping a called ROH
    PLINK.OUT <- NULL
    for(i in unique(plink.combo$id)){                 ## for each individual,
      true.sub <- true.rohs[true.rohs$id == i,]       ## subset true ROH data,
      sub.plink <- plink.combo[plink.combo$id == i,]  ## subset PLINK results
      for(r in 1:nrow(true.sub)){                     ## loop over true ROHs,
        s <- true.sub$start[r]                        ## save true start,
        e <- true.sub$end[r]                          ## and true end.
        ## check for overlaps in PLINK results
        ## called ROHs beginning outside of true ROH, ending inside
        if(nrow(sub.plink[sub.plink$start < s & sub.plink$end >= s & sub.plink$end <= e,]) > 0){
          temp <- sub.plink[sub.plink$start < s & sub.plink$end >= s & sub.plink$end <= e,]
          for(t in 1:nrow(temp)){                                            ## for each called ROH,
            len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
            true.len <- temp$end[t] - s + 1                                  ## length of overlap with the focal true ROH
            save <- c(temp$id[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t], 
                      temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t], 
                      len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t],
                      s, temp$end[t])
            PLINK.OUT <- rbind(PLINK.OUT, save)
          }
        }
        write.table(PLINK.OUT, paste0(rnd,'/',demo,'_PLINK_overlap_results.txt'),
                    sep = '\t', row.names = FALSE, col.names = !file.exists(paste0(rnd,'/',demo,'_PLINK_overlap_results.txt')), 
                    quote = FALSE, append = TRUE)
        PLINK.OUT <- NULL
        ## called ROHs beginning inside of a true ROH, ending outside
        if(nrow(sub.plink[sub.plink$start >= s & sub.plink$start <= e & sub.plink$end > e,]) > 0){
          temp <- sub.plink[sub.plink$start >= s & sub.plink$start <= e & sub.plink$end > e,]
          for(t in 1:nrow(temp)){                                            ## for each called ROH,
            len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
            true.len <- e - temp$start[t] + 1                                ## length of overlap with the focal true ROH
            save <- c(temp$id[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t],
                      temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t],
                      len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t],
                      temp$start[t], e)
            PLINK.OUT <- rbind(PLINK.OUT, save)
          }
        }
        write.table(PLINK.OUT, paste0(rnd,'/',demo,'_PLINK_overlap_results.txt'),
                    sep = '\t', row.names = FALSE, col.names = !file.exists(paste0(rnd,'/',demo,'_PLINK_overlap_results.txt')), 
                    quote = FALSE, append = TRUE)
        PLINK.OUT <- NULL
        ## called ROHs completely covering a true ROH
        if(nrow(sub.plink[sub.plink$start < s & sub.plink$end > e,]) > 0){
          temp <- sub.plink[sub.plink$start < s & sub.plink$end > e,]
          for(t in 1:nrow(temp)){                                            ## for each called ROH,
            len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
            true.len <- e - s + 1                                            ## length of overlap with the focal true ROH,
            save <- c(temp$id[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t],
                      temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t],
                      len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t],
                      s, e)
            PLINK.OUT <- rbind(PLINK.OUT, save)
          }
        }
        write.table(PLINK.OUT, paste0(rnd,'/',demo,'_PLINK_overlap_results.txt'),
                    sep = '\t', row.names = FALSE, col.names = !file.exists(paste0(rnd,'/',demo,'_PLINK_overlap_results.txt')), 
                    quote = FALSE, append = TRUE)
        PLINK.OUT <- NULL
        ## called ROHs completely within a true ROH
        if(nrow(sub.plink[sub.plink$start >= s & sub.plink$end <= e,]) > 0){
          temp <- sub.plink[sub.plink$start >= s & sub.plink$end <= e,]
          for(t in 1:nrow(temp)){                                            ## for each called ROH,
            len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
            true.len <- len                                                  ## length of overlap with the focal true ROH,
            save <- c(temp$id[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t],
                      temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t],
                      len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t],
                      temp$start[t], temp$end[t])
            PLINK.OUT <- rbind(PLINK.OUT, save)
          }
        }
        write.table(PLINK.OUT, paste0(rnd,'/',demo,'_PLINK_overlap_results.txt'),
                    sep = '\t', row.names = FALSE, col.names = !file.exists(paste0(rnd,'/',demo,'_PLINK_overlap_results.txt')), 
                    quote = FALSE, append = TRUE)
        PLINK.OUT <- NULL
      }
    }
    rm(PLINK.OUT)
  }
  
  ##### Read in overlap results for all setting combos #####
  plink.out <- read.table(paste0(rnd,'/',demo,'_PLINK_overlap_results.txt'),
                          sep = '\t', header = TRUE)
  plink.out <- plink.out[plink.out$len >= 100000,]
  
  ## create combo variable to loop over
  plink.out$combo <- paste0(plink.out[,2],'-',
                            plink.out[,3],'-',
                            plink.out[,4],'-',
                            plink.out[,5],'-',
                            plink.out[,6],'-',
                            plink.out[,7],'-',
                            plink.out[,8],'-',
                            plink.out[,9],'-',
                            plink.out[,10])
  
  ##### 3. Calculate true and called f(ROH) values #####
  OUT <- NULL
  for(c in unique(plink.out$combo)){
    sub <- plink.out[plink.out$combo == c,]
      for(i in unique(sub$id)){
        temp <- sub[sub$id == i, c('len','called.roh.id')]
        temp <- temp[!duplicated(temp),]
        true.froh <- sum(true.rohs[true.rohs$id == i, 'length'])/30e6
        call.froh <- sum(temp$len)/30e6
        save <- c(i, unlist(sub[1, c(2:10)]), true.froh, call.froh)
        OUT <- rbind(OUT, save)
      }
    write.table(OUT, paste0(rnd,'/',demo,'_PLINK_individual_froh_results.txt'),
                sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    OUT <- NULL
  }

  OUT <- read.table(paste0(rnd,'/',demo,'_PLINK_individual_froh_results.txt'), header = TRUE)
  froh.stats <- as.data.frame(OUT)
  froh.stats$abs.diff <- abs(froh.stats$call.froh - froh.stats$true.froh)
  froh.stats$combo <- paste0(froh.stats[,2],'-',
                             froh.stats[,3],'-',
                             froh.stats[,4],'-',
                             froh.stats[,5],'-',
                             froh.stats[,6],'-',
                             froh.stats[,7],'-',
                             froh.stats[,8],'-',
                             froh.stats[,9],'-',
                             froh.stats[,10])
  froh.stats <- froh.stats[,c('id','true.froh','call.froh','abs.diff','combo',
                              'covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk')]
  write.csv(froh.stats, paste0(rnd,'/',demo,'_PLINK_individual_froh_results.csv'), row.names = FALSE)
  froh.stats <- read.csv(paste0(rnd,'/',demo,'_PLINK_individual_froh_results.csv'))
  
  ## add f(ROH) difference with sign (i.e., directionality)
  froh.stats$froh.diff <- froh.stats$call.froh - froh.stats$true.froh
  #### if positive, overestimated
  #### if negative, underestimated
  
  
  ##### Writing files for sensitivity analyses (Samarth) #####
  ## calculate population-level data for sensitivity analyses;
  OUT1 <- NULL
  for(g in unique(froh.stats$combo)){
    sub <- froh.stats[froh.stats$combo == g,]
    sub$covg <- gsub('x','', sub$covg)
    if(length(unique(sub$id)) > 40){
      n.indivs <- length(unique(sub$id))
      mean.true <- mean(sub$true.froh)
      mean.called <- mean(sub$call.froh)
      mean.abs.diff <- mean(sub$abs.diff)
      sd.abs.diff <- sd(sub$abs.diff)
      mean.diff <- mean(sub$froh.diff)
      sd.diff <- sd(sub$froh.diff)
      save <- c(n.indivs, mean.true, mean.called, mean.abs.diff, sd.abs.diff, mean.diff, sd.diff, unlist(sub[1,c(5:14)]))
      OUT1 <- rbind(OUT1, save)
    } else{
      next
    }
  }
  colnames(OUT1)[1:7] <- c('n.indivs','mean.true.froh','mean.called.froh','mean.abs.diff','sd.abs.diff','mean.diff','sd.diff')
  pop.lev.res <- as.data.frame(OUT1)
  for(c in c(1:7,9)){
    pop.lev.res[,c] <- as.numeric(pop.lev.res[,c])
  }

  write.csv(OUT1[,c(2,3,7,8:ncol(OUT1))], paste0(rnd,'/',demo,'_population_froh_results_SA.csv'), 
            row.names = FALSE, append = TRUE)
  write.csv(froh.stats[,c(1:3,6:ncol(froh.stats))], paste0(rnd,'/',demo,'_individual_froh_results_SA.csv'), 
            row.names = FALSE, append = TRUE)
}

