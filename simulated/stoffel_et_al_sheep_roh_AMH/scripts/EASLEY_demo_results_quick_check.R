library(scales)

setwd('/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/')

c.len <- 30e6

## bin sizes for length-specific f(ROH) calcs
b1 <- 500e3
b2 <- 1e6
b3 <- 2e6

# dirs <- list.files()
dirs <- c('decline','bottle','small','large-1000','large-2000')

pdf('demo_results_overview.pdf', width = 12, height = 4.5)
par(mfrow = c(1,3))

for(d in dirs){
  print(d)
  setwd(paste0('/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/',d,'/vcfs/'))
  fn <- list.files(pattern = 'final*')
  vcf <- read.table(list.files(pattern = 'final*'))
  colnames(vcf) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', c(paste0('i',c(1:50))))
  
  ##### Check ROH / heterozygosity statistics #####
  ## check ROH length distributions based on variants in VCF
  samp.cols <- grep('i', colnames(vcf)) ## get columns that contain sample data
  
  OUT <- NULL    ## where to save het sites for each individual
  OUT1 <- NULL   ## where to save true ROH coordinates for each individual
  for(col in samp.cols){
    print(colnames(vcf)[col])
    samp <- gsub('i','',colnames(vcf)[col])
    temp <- vcf[,c(2,col)] ## subset to sample genotypes
    temp[temp[,2] == '0|1', 2] <- '1'
    temp[temp[,2] == '1|0', 2] <- '1'
    temp[temp[,2] == '1|1', 2] <- '0'
    temp[temp[,2] == '0|0', 2] <- '0'
    het.sites <- temp[temp[,2] == 1,]
    het.sites <- het.sites[order(het.sites$POS),]
    het.sites[,2] <- samp
    colnames(het.sites) <- c('het.pos','samp')
    OUT <- rbind(OUT, het.sites)
    
    for(h in 1:nrow(het.sites)){
      
      if(nrow(het.sites) == 1){
        s <- 1
        e <- het.sites[1, 1] - 1
        save <- c(samp, s, e, (e-s+1))
        OUT1 <- rbind(OUT1, save)
        
        s <- het.sites[1, 1] + 1
        e <- c.len
        save <- c(samp, s, e, (e-s+1))
        OUT1 <- rbind(OUT1, save)
      } else{
        
        if(h < nrow(het.sites) & het.sites[h, 1] == (het.sites[h+1, 1] + 1)){
          print(paste0(colnames(vcf)[col]," - CONSECUTIVE HET SITES"))
        } else{
          if(h == 1){
            s <- 1
            e <- het.sites[h, 1] - 1
            save <- c(samp, s, e, (e-s+1))
            OUT1 <- rbind(OUT1, save)
          }
          if(h == nrow(het.sites)){
            s <- het.sites[h-1, 1] + 1
            e <- het.sites[h, 1] - 1
            save <- c(samp, s, e, (e-s+1))
            OUT1 <- rbind(OUT1, save)
            s <- het.sites[h, 1] + 1
            e <- c.len
            save <- c(samp, s, e, (e-s+1))
            OUT1 <- rbind(OUT1, save)
          }
          if(h != 1 & h != nrow(het.sites)){
            s <- het.sites[h-1, 1] + 1
            e <- het.sites[h, 1] - 1
            save <- c(samp, s, e, (e-s+1))
            OUT1 <- rbind(OUT1, save)
          }
        }
      }
    }
    if(col == min(samp.cols)){
      write.table(OUT, 'temp_true_het_sites.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
      write.table(OUT1, 'temp_true_roh_coords.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
      OUT <- NULL
      OUT1 <- NULL
    } else{
      write.table(OUT, 'temp_true_het_sites.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
      write.table(OUT1, 'temp_true_roh_coords.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
      OUT <- NULL
      OUT1 <- NULL
    }
  }
  hets <- read.table('temp_true_het_sites.txt')
  colnames(hets) <- c('pos','id')
  
  rohs <- read.table('temp_true_roh_coords.txt')  
  colnames(rohs) <- c('id','s','e','length')
  
  HETS <- NULL
  for(i in unique(hets$id)){
    HETS <- c(HETS, nrow(hets[hets$id == i,])/c.len)
  }
  hist(HETS, main = paste0('Model: ',d,'\nHeterozygosity'), xlab = 'Heterozygosity (# het sites / chrom length)',
       xlim = c(0, 7e-05))
  
  FROH <- NULL
  for(i in unique(rohs$id)){
    sum(rohs[rohs$id == i & rohs$length >= 100e3, 'length'])/c.len
    FROH <- c(FROH, sum(rohs[rohs$id == i & rohs$length >= 100e3, 'length'])/c.len)
  }
  hist(FROH, main = paste0('Model: ',d,'\nf(ROH)'), xlab = 'f(ROH) (ROHs >= 100 kb / chrom length)',
       xlim = c(0, 1))
  
  BINS <- NULL
  for(i in unique(rohs$id)){
    f.bin1 <- sum(rohs[which(rohs$id == i & rohs$length >= 100e3 & rohs$length < b1), 'length'])/c.len
    f.bin2 <- sum(rohs[which(rohs$id == i & rohs$length >= b1 & rohs$length < b2), 'length'])/c.len
    f.bin3 <- sum(rohs[which(rohs$id == i & rohs$length >= b2 & rohs$length < b3), 'length'])/c.len
    f.bin4 <- sum(rohs[which(rohs$id == i & rohs$length >= b3), 'length'])/c.len
    save <- c(i, f.bin1, f.bin2, f.bin3, f.bin4)
    BINS <- rbind(BINS, save)
  }
  plot(0,0, xlim = c(0.75, 4.25), ylim = c(0, 1), col = 'transparent', xaxt = 'n', 
       main = paste0('Model: ',d,'\nBinned f(ROH)'), xlab = '', ylab = 'f(ROH)')
    axis(1, at = c(1:4), labels = c('small','med','large','x-large'))
    points(jitter(rep(1, nrow(BINS)), amount = 0.15), BINS[,2], col = alpha('springgreen4', 0.5), pch = 19)
    points(jitter(rep(2, nrow(BINS)), amount = 0.15), BINS[,3], col = alpha('springgreen4', 0.5), pch = 19)
    points(jitter(rep(3, nrow(BINS)), amount = 0.15), BINS[,4], col = alpha('springgreen4', 0.5), pch = 19)
    points(jitter(rep(4, nrow(BINS)), amount = 0.15), BINS[,5], col = alpha('springgreen4', 0.5), pch = 19)
}
dev.off()
