library(scales)
`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/data/')

##### Set ID numbers of interest #####
f.id <- 'SRR9329986'
m.id <- 'SRR9329990'
o.id <- 'SRR9329977'

##### Loop over ROH-calling results to ID ROHverlap #####
# fns <- list.files(pattern = '*RG_ONLY*')
# for(f in fns){
#   ##### Set up output file
#   ## Output: for each ROH in one parent, identify all overlapping ROHs in the other. Each instance
#   ## of overlap will occupy one line in the output matrix. The output matrix will contain coordinates
#   ## for each parental ROH and the ROHverlap.
#   ROHVERLAP.OUT <- matrix(c('chrom','f.start','f.end','m.start', 'm.end','overlap.start','overlap.end','overlap.len'), nrow = 1)
#   rohverlap.out <- as.data.frame(ROHVERLAP.OUT)
# 
#   print(f)
#   method <- strsplit(f, split = '_')[[1]][2]
#   covg <- strsplit(f, split = '_')[[1]][3]
#   hmm <- strsplit(f, split = '_')[[1]][4]
# 
#   write.table(rohverlap.out, paste0('../output/parental_rohverlap_results_',method,'_',covg,'_',hmm,'.txt'),
#               sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
# 
#   dat <- read.table(f)
#   dat <- dat[,-1]
#   colnames(dat) <- c('id','chrom','start','end','length','n.snps','quality')
#   dat <- dat[dat$length >= 100e3,]
#   f.dat <- dat[dat$id == f.id,]
#   f.dat$roh.id <- c(1:nrow(f.dat))
#   m.dat <- dat[dat$id == m.id,]
#   m.dat$roh.id <- c(1:nrow(m.dat))
# 
#   ## For each paternal ROH, calculate...
#   ## the proportion overlapping a maternal ROH
#   ROHVERLAP.OUT <- NULL
#   for(r in 1:nrow(f.dat)){                     ## loop over paternal ROHs,
#     s <- f.dat$start[r]                        ## save start,
#     e <- f.dat$end[r]                          ## and end.
#     c <- f.dat$chrom[r]
#     ## called ROHs beginning outside of true ROH, ending inside
#     if(nrow(m.dat[m.dat$start < s & m.dat$end >= s & m.dat$end <= e & m.dat$chrom == c,]) > 0){
#       temp <- m.dat[m.dat$start < s & m.dat$end >= s & m.dat$end <= e & m.dat$chrom == c,]
#       for(t in 1:nrow(temp)){                                            ## for each paternal ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- temp$end[t] - s + 1                               ## length of overlap with a maternal ROH
#         save <- c(f.dat$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#     ## called ROHs beginning inside of a true ROH, ending outside
#     if(nrow(m.dat[m.dat$start >= s & m.dat$start <= e & m.dat$end > e & m.dat$chrom == c,]) > 0){
#       temp <- m.dat[m.dat$start >= s & m.dat$start <= e & m.dat$end > e & m.dat$chrom == c,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- e - temp$start[t] + 1                                ## length of overlap with the focal true ROH
#         save <- c(f.dat$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#     ## called ROHs completely covering a true ROH
#     if(nrow(m.dat[m.dat$start < s & m.dat$end > e & m.dat$chrom == c,]) > 0){
#       temp <- m.dat[m.dat$start < s & m.dat$end > e & m.dat$chrom == c,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- e - s + 1                                            ## length of overlap with the focal true ROH,
#         save <- c(f.dat$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#     ## called ROHs completely within a true ROH
#     if(nrow(m.dat[m.dat$start >= s & m.dat$end <= e & m.dat$chrom == c,]) > 0){
#       temp <- m.dat[m.dat$start >= s & m.dat$end <= e & m.dat$chrom == c,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- len                                                  ## length of overlap with the focal true ROH,
#         save <- c(f.dat$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#   }
#   rohverlap.out <- as.data.frame(ROHVERLAP.OUT)
#   write.table(rohverlap.out, paste0('../output/parental_rohverlap_results_',method,'_',covg,'_',hmm,'.txt'),
#               col.names = !file.exists(paste0('../output/parental_rohverlap_results_',method,'_',covg,'_',hmm,'.txt')),
#               sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE)
# 
#   ## plot some random spot-checks to see if script is working
#   colnames(rohverlap.out) <- c('chrom','f.start','f.end','m.start', 'm.end','overlap.start','overlap.end','overlap.len')
#   for(c in c(2:ncol(rohverlap.out))){
#     rohverlap.out[,c] <- as.numeric(rohverlap.out[,c])
#   }
#   samps <- sample(c(1:nrow(rohverlap.out)), 10)
#   w.size <- 5e6
# 
#   pdf(paste0('../figures/parental_rohverlap_results_',method,'_',covg,'_',hmm,'.pdf'), width = 14, height = 5)
#   par(mar = c(5.1, 9.1, 4.1, 2.1))
#   for(s in samps){
#     plot(0, 0, xlim = c(rohverlap.out$overlap.start[s]-(w.size/2), rohverlap.out$overlap.start[s]+(w.size/2)), xlab = paste0('Chromosome position - ',rohverlap.out$chrom[s]),
#          ylim = c(0.7, 3.30), ylab = '', yaxt = 'n')
#       axis(2, at = c(3:1), labels = c('Paternal ROH','Maternal ROH','Parental ROHverlap'), las = 1)
#       sub <- rohverlap.out[rohverlap.out$chrom == rohverlap.out$chrom[s],]
#       y <- 3
#       ## paternal ROHs
#       for(r in 1:nrow(sub)){
#         polygon(c(sub$f.start[r], sub$f.end[r], sub$f.end[r], sub$f.start[r]), c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'dodgerblue2', border = NA)
#       }
#       y <- y-1
#       ## maternal ROHs
#       for(r in 1:nrow(sub)){
#         polygon(c(sub$m.start[r], sub$m.end[r], sub$m.end[r], sub$m.start[r]), c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'dodgerblue2', border = NA)
#       }
#       y <- y-1
#       ## parental ROHverlap
#       for(r in 1:nrow(sub)){
#         lines(c(sub$overlap.start[r], sub$overlap.start[r]), c(y-0.25, 3.25), lty = 3, col = 'darkgrey')
#         lines(c(sub$overlap.end[r], sub$overlap.end[r]), c(y-0.25, 3.25), lty = 3, col = 'darkgrey')
#         polygon(c(sub$overlap.start[r], sub$overlap.end[r], sub$overlap.end[r], sub$overlap.start[r]), c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'dodgerblue4', border = NA)
#       }
#   }
#   dev.off()
# }


##### Loop over results, calculate proportion of parental ROHverlap that's also a ROH in the offspring #####
# setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/output/')
# fns <- list.files(pattern = 'parental_rohverlap_results_*')
# for(f in fns){
#   print(f)
#   method <- strsplit(f, split = '_')[[1]][4]
#   covg <- strsplit(f, split = '_')[[1]][5]
#   hmm <- strsplit(f, split = '_')[[1]][6]
#   hmm <- gsub('.txt', '', hmm)
# 
#   ## read in parental ROHverlap information
#   rohverlap <- read.table(f, header = TRUE)
# 
#   ##### Set up output file
#   ## Output: for each ROHverlap in the parents, identify overlapping ROH(s) in the offspring. Each instance
#   ## of overlap will occupy one line in the output matrix. The output matrix will contain:
#   ## roh.id for paternal ROH, roh.id for maternal ROH, and length of overlap.
#   TRIAD.ROHVERLAP.OUT <- matrix(c('chrom','p.rohv.start','p.rohv.end','o.roh.start', 'o.roh.end','triad.overlap.start','triad.overlap.end','triad.overlap.len'), nrow = 1)
#   triad.rohverlap.out <- as.data.frame(TRIAD.ROHVERLAP.OUT)
# 
#   write.table(triad.rohverlap.out, paste0('../output/parent-offspring_rohverlap_results_',method,'_',covg,'_',hmm,'.txt'),
#               sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
# 
#   ## read in offspring ROH information
#   o.dat <- read.table(paste0('../data/tasdev_',method,'_',covg,'_',hmm,'_RG_ONLY.txt'))
#   o.dat <- o.dat[,-1]
#   colnames(o.dat) <- c('id','chrom','start','end','length','n.snps','quality')
#   o.dat <- o.dat[o.dat$length >= 100e3,]
#   o.dat <- o.dat[o.dat$id == o.id,]
#   o.dat$roh.id <- c(1:nrow(o.dat))
# 
# 
#   ## For each parental ROHverlap segment, calculate...
#   ## the proportion overlapping an offspring ROH
#   TRIAD.ROHVERLAP.OUT <- NULL
#   for(r in 1:nrow(rohverlap)){
#     s <- rohverlap$overlap.start[r]
#     e <- rohverlap$overlap.end[r]
#     c <- rohverlap$chrom[r]
#     ## called ROHs beginning outside of true ROH, ending inside
#     if(nrow(o.dat[o.dat$start < s & o.dat$end >= s & o.dat$end <= e & o.dat$chrom == c,]) > 0){
#       temp <- o.dat[o.dat$start < s & o.dat$end >= s & o.dat$end <= e & o.dat$chrom == c,]
#       for(t in 1:nrow(temp)){
#         len <- temp$end[t] - temp$start[t] + 1
#         overlap.len <- temp$end[t] - s + 1
#         save <- c(rohverlap$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         TRIAD.ROHVERLAP.OUT <- rbind(TRIAD.ROHVERLAP.OUT, save)
#       }
#     }
#     ## called ROHs beginning inside of a true ROH, ending outside
#     if(nrow(o.dat[o.dat$start >= s & o.dat$start <= e & o.dat$end > e & o.dat$chrom == c,]) > 0){
#       temp <- o.dat[o.dat$start >= s & o.dat$start <= e & o.dat$end > e & o.dat$chrom == c,]
#       for(t in 1:nrow(temp)){
#         len <- temp$end[t] - temp$start[t] + 1
#         overlap.len <- e - temp$start[t] + 1
#         save <- c(rohverlap$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         TRIAD.ROHVERLAP.OUT <- rbind(TRIAD.ROHVERLAP.OUT, save)
#       }
#     }
#     ## called ROHs completely covering a true ROH
#     if(nrow(o.dat[o.dat$start < s & o.dat$end > e & o.dat$chrom == c,]) > 0){
#       temp <- o.dat[o.dat$start < s & o.dat$end > e & o.dat$chrom == c,]
#       for(t in 1:nrow(temp)){
#         len <- temp$end[t] - temp$start[t] + 1
#         overlap.len <- e - s + 1
#         save <- c(rohverlap$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         TRIAD.ROHVERLAP.OUT <- rbind(TRIAD.ROHVERLAP.OUT, save)
#       }
#     }
#     o.dat[o.dat$start >= s & o.dat$end <= e,]
#     ## called ROHs completely within a true ROH
#     if(nrow(o.dat[o.dat$start >= s & o.dat$end <= e & o.dat$chrom == c,]) > 0){
#       temp <- o.dat[o.dat$start >= s & o.dat$end <= e & o.dat$chrom == c,]
#       for(t in 1:nrow(temp)){
#         len <- temp$end[t] - temp$start[t] + 1
#         overlap.len <- len
#         save <- c(rohverlap$chrom[r], s, e, temp$start[t], temp$end[t],
#                   max(c(s, temp$start[t])), min(c(e, temp$end[t])), overlap.len)
#         TRIAD.ROHVERLAP.OUT <- rbind(TRIAD.ROHVERLAP.OUT, save)
#       }
#     }
#   }
#   triad.rohverlap.out <- as.data.frame(TRIAD.ROHVERLAP.OUT)
#   write.table(triad.rohverlap.out, paste0('../output/parent-offspring_rohverlap_results_',method,'_',covg,'_',hmm,'.txt'),
#               col.names = !file.exists(paste0('../output/parent-offspring_rohverlap_results_',method,'_',covg,'_',hmm,'.txt')),
#               sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE)
# }

##### Make some random plots to spot-check data summarization #####
# setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/output/')
# fns <- list.files(pattern = 'parent-offspring_rohverlap_results_*')
# for(f in fns){
#   print(f)
#   method <- strsplit(f, split = '_')[[1]][4]
#   covg <- strsplit(f, split = '_')[[1]][5]
#   hmm <- strsplit(f, split = '_')[[1]][6]
#   hmm <- gsub('.txt', '', hmm)
#   
#   ## read in parental andtriad rohverlap results
#   parental.rohverlap <- read.table(paste0('parental_rohverlap_results_',method,'_',covg,'_',hmm,'.txt'))
#   colnames(parental.rohverlap) <- c('chrom','f.start','f.end','m.start', 'm.end','overlap.start','overlap.end','overlap.len')
#   triad.rohverlap <- read.table(f, header = TRUE)
# 
#   ## read in original ROH calls for all 3 individuals
#   dat <- read.table(paste0('../data/tasdev_',method,'_',covg,'_',hmm,'_RG_ONLY.txt'))
#   dat <- dat[,-1]
#   colnames(dat) <- c('id','chrom','start','end','length','n.snps','quality')
#   dat <- dat[dat$length >= 100e3,]
#   # dat <- dat[dat$chrom == 'LR735557.1',]  ## just running samples on 1 chromosome for now to save time (lots of ROHs)
#   f.dat <- dat[dat$id == f.id,]
#   f.dat$roh.id <- c(1:nrow(f.dat))
#   m.dat <- dat[dat$id == m.id,]
#   m.dat$roh.id <- c(1:nrow(m.dat))
#   o.dat <- dat[dat$id == o.id,]
#   o.dat$roh.id <- c(1:nrow(o.dat))
#   
#   ## plot 10 random windows
#   ## set window size
#   w.size <- 5e6
#   
#   pdf(paste0('../figures/triad_rohverlap_results_',method,'_',covg,'_',hmm,'.pdf'), width = 14, height = 5)
#   par(mar = c(5.1, 9.1, 4.1, 2.1))
#   samps <- sample(c(1:nrow(f.dat)), 10)
#   for(s in samps){
#     plot(0, 0, xlim = c(f.dat$start[s]-(w.size/2), f.dat$start[s]+(w.size/2)), xlab = paste0('Chromosome position - ',f.dat$chrom[s]),
#          ylim = c(0.7, 5.30), ylab = '', yaxt = 'n')
#       axis(2, at = c(5:1), labels = c('Paternal ROH','Maternal ROH','Parental ROHverlap','Offspring ROH','Triad ROHverlap'), las = 1)
#       sub <- f.dat[f.dat$chrom == f.dat$chrom[s],]
#       y <- 5
#       ## paternal ROHs
#       for(r in 1:nrow(sub)){
#         polygon(c(sub$start[r], sub$end[r], sub$end[r], sub$start[r]), c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'dodgerblue2', border = NA)
#       }
#       y <- y-1
#       ## maternal ROHs
#       sub <- m.dat[m.dat$chrom == f.dat$chrom[s],]
#       for(r in 1:nrow(sub)){
#         polygon(c(sub$start[r], sub$end[r], sub$end[r], sub$start[r]), c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'dodgerblue2', border = NA)
#       }
#       y <- y-1
#       ## parental ROHverlap
#       sub <- parental.rohverlap[parental.rohverlap$chrom == f.dat$chrom[s],]
#       for(r in 1:nrow(sub)){
#         polygon(c(sub$overlap.start[r], sub$overlap.end[r], sub$overlap.end[r], sub$overlap.start[r]), c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'dodgerblue4', border = NA)
#       }
#       y <- y-1
#       ## offspring ROHs
#       sub <- o.dat[o.dat$chrom == f.dat$chrom[s],]
#       for(r in 1:nrow(sub)){
#         polygon(c(sub$start[r], sub$end[r], sub$end[r], sub$start[r]), c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'goldenrod1', border = NA)
#       }
#       ## triad ROHverlap
#       y <- y-1
#       sub <- triad.rohverlap[triad.rohverlap$chrom == f.dat$chrom[s],]
#       for(r in 1:nrow(sub)){
#         polygon(c(sub$triad.overlap.start[r], sub$triad.overlap.end[r], sub$triad.overlap.end[r], sub$triad.overlap.start[r]), 
#                 c(y-0.25, y-0.25, y+0.25, y+0.25), lwd = 4, col = 'springgreen3', border = NA)
#       }
#   }
#   dev.off()
# }


##### Calculate proportions of parental ROHverlap called as ROH in offspring #####
setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/output/')

CALC.PROPS <- matrix(c('method','covg','hmm','p.roh.id','triad.roh.id','p.rohv.len','triad.rohv.len','prop.p.cov'), nrow = 1)

fns <- list.files(pattern = 'parent-offspring_rohverlap_results_*')
for(f in fns){
  print(f)
  method <- strsplit(f, split = '_')[[1]][4]
  covg <- strsplit(f, split = '_')[[1]][5]
  hmm <- strsplit(f, split = '_')[[1]][6]
  hmm <- gsub('.txt', '', hmm)

  ## read in parental and triad rohverlap results
  parental.rohverlap <- read.table(paste0('parental_rohverlap_results_',method,'_',covg,'_',hmm,'.txt'), header = TRUE)
  parental.rohverlap$p.roh.id <- c(1:nrow(parental.rohverlap))
  triad.rohverlap <- read.table(f, header = TRUE)
  triad.rohverlap$t.roh.id <- c(1:nrow(triad.rohverlap))

  ## read in original ROH calls for all 3 individuals
  dat <- read.table(paste0('../data/tasdev_',method,'_',covg,'_',hmm,'_RG_ONLY.txt'))
  dat <- dat[,-1]
  colnames(dat) <- c('id','chrom','start','end','length','n.snps','quality')
  dat <- dat[dat$length >= 100e3,]
  f.dat <- dat[dat$id == f.id,]
  f.dat$roh.id <- c(1:nrow(f.dat))
  m.dat <- dat[dat$id == m.id,]
  m.dat$roh.id <- c(1:nrow(m.dat))
  o.dat <- dat[dat$id == o.id,]
  o.dat$roh.id <- c(1:nrow(o.dat))
  
  ## loop over parental ROHverlap results, checking to see if the overlap was caught in the triad overlap results.
  ## if it was, calculate the length of overlap between the two sets of results.
  ## if it wasn't, record its length as a miss. 
  ## for now, don't filter based on overlap length, that can be done downstream.
  for(r in 1:nrow(parental.rohverlap)){
    p.c <- parental.rohverlap$chrom[r]
    p.s <- parental.rohverlap$overlap.start[r]
    p.e <- parental.rohverlap$overlap.end[r]
    p.roh.id <- parental.rohverlap$p.roh.id[r]
    
    ## see if it exists in the triad rohverlap results, exactly once
    if(nrow(triad.rohverlap[triad.rohverlap$chrom == p.c &
                            triad.rohverlap$p.rohv.start == p.s &
                            triad.rohverlap$p.rohv.end == p.e,]) == 1){
      ## if it does, save and do some calcs
      temp <- triad.rohverlap[triad.rohverlap$chrom == p.c &
                              triad.rohverlap$p.rohv.start == p.s &
                              triad.rohverlap$p.rohv.end == p.e,]
      p.rohv.len <- temp$p.rohv.end - temp$p.rohv.start + 1
      save <- c(method, covg, hmm, p.roh.id, temp$t.roh.id[1], 
                p.rohv.len, temp$triad.overlap.len[1], (p.rohv.len - temp$triad.overlap.len[1])/p.rohv.len)
      CALC.PROPS <- rbind(CALC.PROPS, save)
      ## if it doesn't, do some calcs
    }
    ## see if it exists in the triad rohverlap results, more than once
    if(nrow(triad.rohverlap[triad.rohverlap$chrom == p.c &
                            triad.rohverlap$p.rohv.start == p.s &
                            triad.rohverlap$p.rohv.end == p.e,]) > 1){
      ## if it does, save and do some calcs
      temp <- triad.rohverlap[triad.rohverlap$chrom == p.c &
                                triad.rohverlap$p.rohv.start == p.s &
                                triad.rohverlap$p.rohv.end == p.e,]
      for(t in 1:nrow(temp)){
        p.rohv.len <- temp$p.rohv.end[t] - temp$p.rohv.start[t] + 1
        save <- c(method, covg, hmm, p.roh.id, temp$t.roh.id[t], 
                  p.rohv.len, temp$triad.overlap.len[t], (p.rohv.len - temp$triad.overlap.len[t])/p.rohv.len)
        CALC.PROPS <- rbind(CALC.PROPS, save)
      }
    ## if it doesn't exist in the triad rohvelap results, do some calcs
    } 
    if(nrow(triad.rohverlap[triad.rohverlap$chrom == p.c &
                            triad.rohverlap$p.rohv.start == p.s &
                            triad.rohverlap$p.rohv.end == p.e,]) == 0){
      p.rohv.len <- parental.rohverlap$overlap.len[r]
      triad.rohv.len <- 0
      save <- c(method, covg, hmm, p.roh.id, NA, p.rohv.len, triad.rohv.len, (p.rohv.len - triad.rohv.len)/p.rohv.len)
      CALC.PROPS <- rbind(CALC.PROPS, save)
    }
  }
}
prop.calcs <- as.data.frame(CALC.PROPS)
colnames(prop.calcs) <- prop.calcs[1,]
prop.calcs <- prop.calcs[-1,]
rm(CALC.PROPS)

prop.calcs$combo <- paste0(prop.calcs$method,'-',prop.calcs$covg,'-',prop.calcs$hmm)
for(c in c(6:7)){
  prop.calcs[,c] <- as.numeric(prop.calcs[,c])
}

pdf('../figures/parental_vs_triad_rohverlap_lengths.pdf', width = 10, height = 5)
par(mfrow = c(1,2))
for(c in unique(prop.calcs$combo)){
  sub <- prop.calcs[prop.calcs$combo == c,]
  plot(sub$p.rohv.len, sub$triad.rohv.len, pch = 19, col = alpha('springgreen3', 0.5), 
       main = c, xlab = 'Parental ROHverlap length', ylab = 'Triad ROHverlap length')
    abline(0, 1, lty = 2)
}
dev.off()

OUT <- NULL
for(c in unique(prop.calcs$combo)){
  sub <- prop.calcs[prop.calcs$combo == c,]
  sub.filt <- sub[sub$p.rohv.len >= 100e3,]
  save <- c(c, sum(sub$p.rohv.len), sum(sub$triad.rohv.len), sum(sub$triad.rohv.len)/sum(sub$p.rohv.len),
            sum(sub.filt$p.rohv.len), sum(sub.filt$triad.rohv.len), sum(sub.filt$triad.rohv.len)/sum(sub.filt$p.rohv.len))
  OUT <- rbind(OUT, save)
}
stats <- as.data.frame(OUT)
for(c in c(2:7)){
  stats[,c] <- as.numeric(stats[,c])
}
colnames(stats) <- c('combo','p.rohv.len','triad.rohv.len','prop.p.cov','filt.p.roh.len','filt.triad.rohv.len','prop.filt.p.cov')
stats <- stats[order(stats$prop.p.cov),]
