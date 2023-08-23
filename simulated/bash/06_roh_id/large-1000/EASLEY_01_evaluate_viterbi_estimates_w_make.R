`%notin%` <- Negate(`%in%`)

setwd('/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/data/06_roh_id/output/')

## filename format: large-1000_sample_cvg_${c}_PL_VT_ONLY.txt
fns <- list.files()
fns <- fns[grep('large-1000', fns)]
fns <- fns[grep('VT_ONLY', fns)]

## get coverages and viterbi-training values tested
covgs <- unique(do.call(rbind, strsplit(fns, split = '_'))[,4])

for(d in c('GT','PL')){
  OUT <- NULL
  for(c in covgs){
    dat <- read.table(paste0('large-1000_sample_cvg_',c,'_',d,'_VT_ONLY.txt'))
    dat <- dat[,-c(1,6,8)]
    colnames(dat) <- c('sample','iteration','dAZ','dHW','P(AZ|HW)','P(HW|AZ)')
    for(i in unique(dat$sample)){
      sub <- dat[dat$sample == i,]
      save <- c(i, c, sub$`P(AZ|HW)`[nrow(sub)], sub$`P(HW|AZ)`[nrow(sub)])
      OUT <- rbind(OUT, save)
    }
  }
  OUT1 <- NULL
  for(c in unique(OUT[,2])){
      sub <- OUT[OUT[,2] == c,]
      save <- c(c, mean(as.numeric(sub[,3])), sd(as.numeric(sub[,3])), mean(as.numeric(sub[,4])), sd(as.numeric(sub[,4])))
      OUT1 <- rbind(OUT1, save)
  }
  write.table(OUT1, paste0('large-1000_',d,'_viterbi_estimates.txt'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}