setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/')

c.len <- 30e6

dirs <- list.files()

pdf('demo_results_overview.pdf', width = 10, height = 5)
par(mfrow = c(1,2))

for(d in dirs){
  print(d)
  setwd(paste0('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/',d,'/vcfs/'))
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
  hist(HETS, main = paste0('Model: ',d,'\nHeterozygosity'), xlab = 'Heterozygosity (# het sites / chrom length)')
  
  FROH <- NULL
  for(i in unique(rohs$id)){
    sum(rohs[rohs$id == i & rohs$length >= 100e3, 'length'])/c.len
    FROH <- c(FROH, sum(rohs[rohs$id == i & rohs$length >= 100e3, 'length'])/c.len)
  }
  hist(FROH, main = paste0('Model: ',d,'\nf(ROH)'), xlab = 'f(ROH) (ROHs >= 100 kb / chrom length)',
       xlim = c(0, 1))
}
dev.off()
