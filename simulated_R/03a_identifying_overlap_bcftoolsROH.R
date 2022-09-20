`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/')

##### Read in known/true heterozygosity + ROH information #####
true.rohs <- read.table('true_indiv_data/true_roh_coords.txt')
colnames(true.rohs) <- c('id','start','end','length')
true.rohs <- true.rohs[true.rohs$length >= 100000,] ### only keeping true ROHs >= 100 kb because that's all we're evaluating the ability to call
## create unique ROH ID for linking true ROHs to called ROHs
ns <- c(1:nrow(true.rohs))
true.rohs$true.roh.id <- ns
chrom.len <- 30e6

## from ROH data, get heterozygous sites
true.het <- read.table('true_indiv_data/true_het_sites.txt')
colnames(true.het) <- c('pos','id')

##### Read in ROH calling results #####
## BCFtools/ROH (written from clean_plotting_bcftoolsROH_output.R)
## GT
bcf.gt.res <- read.table('bcftoolsroh_results/bcftoolsROH_GT_all_coordinates.txt',
                         header = TRUE, sep = '\t')
bcf.gt.res <- bcf.gt.res[bcf.gt.res$length >= 100000,] ## applying 100kb filter, as would normally be done
bcf.gt.res$called.roh.id <- c(1:nrow(bcf.gt.res))
## PL
bcf.pl.res <- read.table('bcftoolsroh_results/bcftoolsROH_PL_all_coordinates.txt',
                         header=TRUE, sep = '\t')
bcf.pl.res <- bcf.pl.res[bcf.pl.res$length >= 100000,] ## applying 100kb filter, as would normally be done
bcf.pl.res$called.roh.id <- c(1:nrow(bcf.pl.res))

##### 2A. ROH identification -- BCFtools/ROH results #####
## Output: for each true ROH, identify all overlapping called ROHs. Each instance
## of overlap will occupy one line in the output matrix. The output matrix will contain:
## roh.id for true ROH, id/pop.size/covg/start/end for called ROH, and length of overlap.

## For each true ROH, calculate...
## the proportion overlapping a called ROH
BCF.GT.OUT <- NULL
BCF.PL.OUT <- NULL

for(i in unique(true.rohs$id)){              ## for each individual,
  print(i)
  true.sub <- true.rohs[true.rohs$id == i,]  ## subset true ROH data,
  sub.gt <- bcf.gt.res[bcf.gt.res$id == i,]  ## subset GT results
  sub.pl <- bcf.pl.res[bcf.pl.res$id == i,]  ## subset PL results
  for(r in 1:nrow(true.sub)){                     ## loop over true ROHs,
    s <- true.sub$start[r]                        ## save true start,
    e <- true.sub$end[r]                          ## and true end.
    ## check for overlaps in both sets of BCFtools results
    ## GT
    ## called ROHs beginning outside of true ROH, ending inside
    if(nrow(sub.gt[sub.gt$start < s & sub.gt$end >= s & sub.gt$end <= e,]) > 0){
      temp <- sub.gt[sub.gt$start < s & sub.gt$end >= s & sub.gt$end <= e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- temp$end[t] - s + 1                                  ## length of overlap with the focal true ROH
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.GT.OUT <- rbind(BCF.GT.OUT, save)
      }
    }
    ## called ROHs beginning inside of a true ROH, ending outside
    if(nrow(sub.gt[sub.gt$start >= s & sub.gt$start <= e & sub.gt$end > e,]) > 0){
      temp <- sub.gt[sub.gt$start >= s & sub.gt$start <= e & sub.gt$end > e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- e - temp$start[t] + 1                                ## length of overlap with the focal true ROH
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.GT.OUT <- rbind(BCF.GT.OUT, save)
      }
    }
    ## called ROHs completely covering a true ROH
    if(nrow(sub.gt[sub.gt$start < s & sub.gt$end > e,]) > 0){
      temp <- sub.gt[sub.gt$start < s & sub.gt$end > e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- e - s + 1                                            ## length of overlap with the focal true ROH,
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.GT.OUT <- rbind(BCF.GT.OUT, save)
      }
    }
    sub.gt[sub.gt$start >= s & sub.gt$end <= e,]
    ## called ROHs completely within a true ROH
    if(nrow(sub.gt[sub.gt$start >= s & sub.gt$end <= e,]) > 0){
      temp <- sub.gt[sub.gt$start >= s & sub.gt$end <= e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- len                                                  ## length of overlap with the focal true ROH,
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.GT.OUT <- rbind(BCF.GT.OUT, save)
      }
    }
    ## PL
    ## called ROHs beginning outside of true ROH, ending inside
    if(nrow(sub.pl[sub.pl$start < s & sub.pl$end >= s & sub.pl$end <= e,]) > 0){
      temp <- sub.pl[sub.pl$start < s & sub.pl$end >= s & sub.pl$end <= e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- temp$end[t] - s + 1                                  ## length of overlap with the focal true ROH
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.PL.OUT <- rbind(BCF.PL.OUT, save)
      }
    }
    ## called ROHs beginning inside of a true ROH, ending outside
    if(nrow(sub.pl[sub.pl$start >= s & sub.pl$start <= e & sub.pl$end > e,]) > 0){
      temp <- sub.pl[sub.pl$start >= s & sub.pl$start <= e & sub.pl$end > e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- e - temp$start[t] + 1                                ## length of overlap with the focal true ROH
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.PL.OUT <- rbind(BCF.PL.OUT, save)
      }
    }
    ## called ROHs completely covering a true ROH
    if(nrow(sub.pl[sub.pl$start < s & sub.pl$end > e,]) > 0){
      temp <- sub.pl[sub.pl$start < s & sub.pl$end > e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- e - s + 1                                            ## length of overlap with the focal true ROH,
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.PL.OUT <- rbind(BCF.PL.OUT, save)
      }
    }
    ## called ROHs completely within a true ROH
    if(nrow(sub.pl[sub.pl$start >= s & sub.pl$end <= e,]) > 0){
      temp <- sub.pl[sub.pl$start >= s & sub.pl$end <= e,]
      for(t in 1:nrow(temp)){                                            ## for each called ROH,
        len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
        true.len <- len                                                  ## length of overlap with the focal true ROH,
        save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t],
                             temp$num.snps[t], temp$quality[t]))
        BCF.PL.OUT <- rbind(BCF.PL.OUT, save)
      }
    }
  }
}

colnames(BCF.PL.OUT) <- c('id','pop.size','covg','len','true.len','true.roh.id','called.roh.id','num.called.snps','roh.call.quality')
colnames(BCF.GT.OUT) <- c('id','pop.size','covg','len','true.len','true.roh.id','called.roh.id','num.called.snps','roh.call.quality')
pl.out <- as.data.frame(BCF.PL.OUT)
gt.out <- as.data.frame(BCF.GT.OUT)
write.table(pl.out, 'bcftoolsroh_results/bcftoolsROH_PL_overlap_results.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(gt.out, 'bcftoolsroh_results/bcftoolsROH_GT_overlap_results.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)