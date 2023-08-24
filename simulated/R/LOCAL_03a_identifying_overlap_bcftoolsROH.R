`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/')

##### Read in all ROH-calling results #####
res <- read.table('bcftools_output/bcftoolsROH_all_coordinates.txt', header = TRUE)

##### Loop over demographic scenarios #####
true.fns <- list.files(path = 'slim_true_data/', pattern = 'roh_coords')

for(t in true.fns){
  demo <- unlist(strsplit(t, split = '_'))[1]
  demo.res <- res[res$demo == demo,]

  ##### Read in known/true heterozygosity + ROH information #####
  true.rohs <- read.table(paste0('slim_true_data/',t))
  colnames(true.rohs) <- c('id','start','end','length','true.roh.id')
  true.rohs <- true.rohs[true.rohs$length >= 100000,] ### only keeping true ROHs >= 100 kb because that's all we're evaluating the ability to call
  chrom.len <- 30e6
  
  BCF.OUT <- matrix(c('method','covg','hmm','id','len','true.len','true.roh.id',
                      'called.roh.id','num.called.snps','roh.call.quality'), nrow = 1)
  
  ##### Loop over method / covg / HMM settings within scenario #####
  for(m in unique(demo.res$method)){
    for(c in unique(demo.res$covg)){
      for(h in unique(demo.res$hmm)){
        sub.demo <- demo.res[demo.res$method == m & demo.res$covg == c & demo.res$hmm == h,]
        sub.demo$id <- gsub('i','', sub.demo$id)
        
        ##### Identify overlap between true and called ROHs #####
        ## Output: for each true ROH, identify all overlapping called ROHs. Each instance
        ## of overlap will occupy one line in the output matrix. The output matrix will contain:
        ## roh.id for true ROH, id/pop.size/covg/start/end for called ROH, and length of overlap.
        
        for(i in unique(true.rohs$id)){              ## for each individual,
          print(i)
          true.sub <- true.rohs[true.rohs$id == i,]       ## subset true ROH data,
          sub.bcf <- sub.demo[sub.demo$id == i,]          ## subset BCFtools results
          for(r in 1:nrow(true.sub)){                     ## loop over true ROHs,
            s <- true.sub$start[r]                        ## save true start,
            e <- true.sub$end[r]                          ## and true end.
            ## called ROHs beginning outside of true ROH, ending inside
            if(nrow(sub.bcf[sub.bcf$start < s & sub.bcf$end >= s & sub.bcf$end <= e,]) > 0){
              temp <- sub.bcf[sub.bcf$start < s & sub.bcf$end >= s & sub.bcf$end <= e,]
              for(t in 1:nrow(temp)){                                            ## for each called ROH,
                len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
                true.len <- temp$end[t] - s + 1                                  ## length of overlap with the focal true ROH
                save <- c(m, c, h, temp$id[t], len, true.len, true.sub$true.roh.id[r], 
                                     temp$called.roh.id[t], temp$num.snps[t], temp$quality[t])
                BCF.OUT <- rbind(BCF.OUT, save)
              }
            }
            ## called ROHs beginning inside of a true ROH, ending outside
            if(nrow(sub.bcf[sub.bcf$start >= s & sub.bcf$start <= e & sub.bcf$end > e,]) > 0){
              temp <- sub.bcf[sub.bcf$start >= s & sub.bcf$start <= e & sub.bcf$end > e,]
              for(t in 1:nrow(temp)){                                            ## for each called ROH,
                len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
                true.len <- e - temp$start[t] + 1                                ## length of overlap with the focal true ROH
                save <- c(m, c, h, temp$id[t], len, true.len, true.sub$true.roh.id[r], 
                                     temp$called.roh.id[t], temp$num.snps[t], temp$quality[t])
                BCF.OUT <- rbind(BCF.OUT, save)
              }
            }
            ## called ROHs completely covering a true ROH
            if(nrow(sub.bcf[sub.bcf$start < s & sub.bcf$end > e,]) > 0){
              temp <- sub.bcf[sub.bcf$start < s & sub.bcf$end > e,]
              for(t in 1:nrow(temp)){                                            ## for each called ROH,
                len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
                true.len <- e - s + 1                                            ## length of overlap with the focal true ROH,
                save <- c(m, c, h, temp$id[t], len, true.len, true.sub$true.roh.id[r], 
                                     temp$called.roh.id[t], temp$num.snps[t], temp$quality[t])
                BCF.OUT <- rbind(BCF.OUT, save)
              }
            }
            ## called ROHs completely within a true ROH
            if(nrow(sub.bcf[sub.bcf$start >= s & sub.bcf$end <= e,]) > 0){
              temp <- sub.bcf[sub.bcf$start >= s & sub.bcf$end <= e,]
              for(t in 1:nrow(temp)){                                            ## for each called ROH,
                len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
                true.len <- len                                                  ## length of overlap with the focal true ROH,
                save <- c(m, c, h, temp$id[t], len, true.len, true.sub$true.roh.id[r], 
                                     temp$called.roh.id[t], temp$num.snps[t], temp$quality[t])
                BCF.OUT <- rbind(BCF.OUT, save)
              }
            }
          }
        }
      }
    }
  }
  bcf.out <- as.data.frame(BCF.OUT)
  write.table(bcf.out, paste0('bcftools_output/bcftoolsROH_',demo,'_overlap_results.txt'),
              sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
}
