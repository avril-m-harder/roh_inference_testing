`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/data/')

##### Set ID numbers of interest #####
f.id <- 'SRR9329986'
m.id <- 'SRR9329990'
o.id <- 'SRR9329977'

##### Set up output file #####
## Output: for each ROH in one parent, identify all overlapping ROHs in the other. Each instance
## of overlap will occupy one line in the output matrix. The output matrix will contain:
## roh.id for paternal ROH, roh.id for maternal ROH, and length of overlap.
# ROHVERLAP.OUT <- matrix(c('method','covg','hmm','chrom','f.start','f.end','m.start', 'm.end','overlap.len'), nrow = 1)
# rohverlap.out <- as.data.frame(ROHVERLAP.OUT)
# write.table(rohverlap.out, '../output/triad_rohverlap_results.txt',
#             sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
# 
# ##### Loop over ROH-calling results to ID ROHverlap #####
# fns <- list.files(pattern = '*RG_ONLY*')
# for(f in fns){
#   print(f)
#   method <- strsplit(f, split = '_')[[1]][2]
#   covg <- strsplit(f, split = '_')[[1]][3]
#   hmm <- strsplit(f, split = '_')[[1]][4]
#   
#   dat <- read.table(f)
#   dat <- dat[,-1]
#   colnames(dat) <- c('id','chrom','start','end','length','n.snps','quality')
#   dat <- dat[dat$length >= 100e3,]
#   # dat <- dat[dat$chrom == 'LR735557.1',]  ## just running samples on 1 chromosome for now to save time (lots of ROHs)
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
#     ## called ROHs beginning outside of true ROH, ending inside
#     if(nrow(m.dat[m.dat$start < s & m.dat$end >= s & m.dat$end <= e,]) > 0){
#       temp <- m.dat[m.dat$start < s & m.dat$end >= s & m.dat$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each paternal ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- temp$end[t] - s + 1                               ## length of overlap with a maternal ROH
#         save <- c(method, covg, hmm, f.dat$chrom[r], s, e, temp$start[t], temp$end[t], overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#     ## called ROHs beginning inside of a true ROH, ending outside
#     if(nrow(m.dat[m.dat$start >= s & m.dat$start <= e & m.dat$end > e,]) > 0){
#       temp <- m.dat[m.dat$start >= s & m.dat$start <= e & m.dat$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- e - temp$start[t] + 1                                ## length of overlap with the focal true ROH
#         save <- c(method, covg, hmm, f.dat$chrom[r], s, e, temp$start[t], temp$end[t], overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#     ## called ROHs completely covering a true ROH
#     if(nrow(m.dat[m.dat$start < s & m.dat$end > e,]) > 0){
#       temp <- m.dat[m.dat$start < s & m.dat$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- e - s + 1                                            ## length of overlap with the focal true ROH,
#         save <- c(method, covg, hmm, f.dat$chrom[r], s, e, temp$start[t], temp$end[t], overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#     m.dat[m.dat$start >= s & m.dat$end <= e,]
#     ## called ROHs completely within a true ROH
#     if(nrow(m.dat[m.dat$start >= s & m.dat$end <= e,]) > 0){
#       temp <- m.dat[m.dat$start >= s & m.dat$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         overlap.len <- len                                                  ## length of overlap with the focal true ROH,
#         save <- c(method, covg, hmm, f.dat$chrom[r], s, e, temp$start[t], temp$end[t], overlap.len)
#         ROHVERLAP.OUT <- rbind(ROHVERLAP.OUT, save)
#       }
#     }
#   rohverlap.out <- as.data.frame(ROHVERLAP.OUT)
#   rohverlap.out <- rohverlap.out[!duplicated(rohverlap.out),]
#   write.table(rohverlap.out, '../output/triad_rohverlap_results.txt',
#               col.names = !file.exists('../output/triad_rohverlap_results.txt'),
#               sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE)
#   }
# }

rohverlap <- read.table('../output/triad_rohverlap_results.txt', header = TRUE)
rohverlap <- rohverlap[!duplicated(rohverlap),]

##### Randomly select a few instances of ROHverlap and plot data for all 3 individuals #####
samps <- sample(c(1:nrow(rohverlap)), 10)

for(s in samps){
  to.samp <- rohverlap[s,]
  dat <- read.table(paste0('tasdev_',to.samp$method[1],'_',to.samp$covg[1],'_',to.samp$hmm[1],'_RG_ONLY.txt'))
  dat <- dat[,-1]
  colnames(dat) <- c('id','chrom','start','end','length','n.snps','quality')
  dat <- dat[dat$length >= 100e3,]
  dat <- dat[dat$chrom == to.samp$chrom[1],]
  f.dat <- dat[dat$id == f.id,]
  m.dat <- dat[dat$id == m.id,]
  o.dat <- dat[dat$id == o.id,]
  
  plot(0, 0, xlim = c(min(to.samp[,c(5:8)])-100e3, max(to.samp[,c(5:8)])+100e3), ylim = c(1,3),
       xlab = 'Chromosome position', ylab = '', yaxt = 'n',
       main = paste0(to.samp$method[1],' - ',to.samp$covg[1],' - ',to.samp$hmm[1]))
    axis(2, at = c(1:3), labels = c('Offspring','Mother','Father'))
    for(r in 1:nrow(f.dat)){
      lines(c(f.dat$start[r], f.dat$end[r]), c(3,3), lwd = 6, col = 'dodgerblue3')
    }
    for(r in 1:nrow(m.dat)){
      lines(c(m.dat$start[r], m.dat$end[r]), c(2,2), lwd = 6, col = 'dodgerblue3')
    }
    for(r in 1:nrow(o.dat)){
      lines(c(o.dat$start[r], o.dat$end[r]), c(1,1), lwd = 6, col = 'dodgerblue3')
    }
    
}
