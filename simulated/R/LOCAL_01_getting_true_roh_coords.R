setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/slim_true_data/')
library(ape)
library(scales)
`%notin%` <- Negate(`%in%`)

##### Read in SLiM results #####
## set chromosome length for SLiM run
c.len <- 30e6

fns <- list.files()
## for each demograhpic scenario,
## read in VCF
for(f in fns){
  vcf <- read.table(f, header=TRUE)
  demo <- strsplit(f, split = '_')[[1]][2]
  head(vcf[,c(1:10)])
  
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
    if(col == min(samp.cols)){
      write.table(OUT, paste0(demo,'_true_het_sites.txt'), sep = '\t', row.names = FALSE, 
                  quote = FALSE, col.names = FALSE)
      write.table(OUT1, paste0(demo,'_true_roh_coords.txt'), sep = '\t', row.names = FALSE, 
                  quote = FALSE, col.names = FALSE)
      OUT <- NULL
      OUT1 <- NULL
    } else{
      write.table(OUT, paste0(demo,'_true_het_sites.txt'), sep = '\t', row.names = FALSE, 
                  quote = FALSE, col.names = FALSE, append = TRUE)
      write.table(OUT1, paste0(demo,'_true_roh_coords.txt'), sep = '\t', row.names = FALSE, 
                  quote = FALSE, col.names = FALSE, append = TRUE)
      OUT <- NULL
      OUT1 <- NULL
    }
  }
}


