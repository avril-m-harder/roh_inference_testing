setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/slim_true_data/')
library(ape)
library(scales)
`%notin%` <- Negate(`%in%`)

##### Read in SLiM results #####
## set chromosome length for SLiM run
c.len <- 30e6

fns <- list.files(pattern = '.vcf')
## for each demograhpic scenario,
## read in VCF
for(f in fns){
  vcf <- read.table(f, header=TRUE)
  demo <- strsplit(f, split = '_')[[1]][2]
  print(demo)
  head(vcf[,c(1:10)])
  
  true.roh.id <- 1
  
  ##### Check ROH / heterozygosity statistics #####
  ## check ROH length distributions based on variants in VCF
  samp.cols <- grep('i', colnames(vcf)) ## get columns that contain sample data
  
  OUT <- NULL    ## where to save het sites for each individual
  OUT1 <- NULL   ## where to save true ROH coordinates for each individual
  for(col in samp.cols){
    # print(colnames(vcf)[col])
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
        if(h == 1 & nrow(het.sites)!= 1){
          s <- 1
          e <- het.sites[h, 1] - 1
          save <- c(samp, s, e, (e-s+1), true.roh.id)
          true.roh.id <- true.roh.id + 1
          OUT1 <- rbind(OUT1, save)
        }
        if(h == nrow(het.sites) & nrow(het.sites)!= 1){
          s <- het.sites[h-1, 1] + 1
          e <- het.sites[h, 1] - 1
          save <- c(samp, s, e, (e-s+1), true.roh.id)
          true.roh.id <- true.roh.id + 1
          OUT1 <- rbind(OUT1, save)
          s <- het.sites[h, 1] + 1
          e <- c.len
          save <- c(samp, s, e, (e-s+1), true.roh.id)
          true.roh.id <- true.roh.id + 1
          OUT1 <- rbind(OUT1, save)
        }
        if(h != 1 & h != nrow(het.sites) & nrow(het.sites)!= 1){
          s <- het.sites[h-1, 1] + 1
          e <- het.sites[h, 1] - 1
          save <- c(samp, s, e, (e-s+1), true.roh.id)
          true.roh.id <- true.roh.id + 1
          OUT1 <- rbind(OUT1, save)
        }
        if(nrow(het.sites) == 1){
          s <- 1
          e <- het.sites[h, 1] - 1
          save <- c(samp, s, e, (e-s+1), true.roh.id)
          true.roh.id <- true.roh.id + 1
          OUT1 <- rbind(OUT1, save)
          
          s <- het.sites[h, 1] + 1
          e <- c.len
          save <- c(samp, s, e, (e-s+1), true.roh.id)
          true.roh.id <- true.roh.id + 1
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

### Format histograms for manuscript figures
pdf('/Users/Avril/Desktop/scenario_froh.pdf', width = 5, height = 5)

## histograms require demo-specific # breaks for uniformity
breaks <- c(10, 11, 8, 7)
demos <- c('bottle','decline','small','large-1000')
cex.axis <- 1.5
lwd <- 2
line<-par(lwd = 2)

for(d in 1:4){
  dat <- read.table(paste0(demos[d],'_true_roh_coords.txt'), header = FALSE)
  colnames(dat) <- c('id','start','end','length')
  dat <- dat[dat$length >= 100e3,]
  SAVE <- NULL
  for(i in unique(dat$id)){
    froh <- sum(dat[dat$id == i, 'length'])/c.len
    SAVE <- c(SAVE, froh)
  }
  hist(SAVE, xlim = c(0,1), xlab = '', main = demos[d], ylim = c(0,25), breaks = breaks[d], xaxt = 'n',
       lwd = lwd, ylab = '', yaxt = 'n') ## 19
    axis(1, at = c(0, .25, .5, .75, 1), labels = c('0.0','','0.5','','1.0'), lwd = lwd, cex.axis = cex.axis)
    # axis(2, at = c(0, 5, 10, 15, 20, 25), labels = c('0','','10','','20',''), lwd = 2)
    axis(2, at = c(0, 5, 10, 15, 20, 25), lwd = lwd, cex.axis = cex.axis)
    
  plot(density(SAVE), xlim = c(0, 1), ylim = c(0, 10), col = 'springgreen4', main = demos[d], xaxt = 'n', yaxt = 'n', bty = 'n', zero.line = FALSE)
    axis(1, at = c(0, .25, .5, .75, 1), labels = c('0.0','','0.5','','1.0'), lwd = lwd, cex.axis = cex.axis)
    axis(2, at = c(0, 5, 10), lwd = lwd, cex.axis = cex.axis)
    polygon(density(SAVE), col = 'springgreen4', border = NA)
}

dev.off()


## Length bin figures
b.1 <- 100e3
b.2 <- 500e3
b.3 <- 1e6
b.4 <- 2e6

pdf('/Users/Avril/Desktop/scenario_froh_length_bins.pdf', width = 5, height = 5)
pt.cex <- 1.25
pt.alph <- 0.3
line.len <- 0.3
cex.axis <- 1.5
lwd <- 2

for(f in fns){
  demo <- strsplit(f, split = '_')[[1]][2]
  
  dat <- read.table(paste0(demo,'_true_roh_coords.txt'), header = FALSE)
  colnames(dat) <- c('id','start','end','length')
  dat <- dat[dat$length >= 100e3,]
  OUT <- NULL
  for(i in unique(dat$id)){
    sub <- dat[dat$id == i,]
    b1 <- sum(sub[sub$length >= b.1 & sub$length < b.2, 'length'])/c.len
    b2 <- sum(sub[sub$length >= b.2 & sub$length < b.3, 'length'])/c.len
    b3 <- sum(sub[sub$length >= b.3 & sub$length < b.4, 'length'])/c.len
    b4 <- sum(sub[sub$length >= b.4, 'length'])/c.len
    
    save <- c(i, b1, b2, b3, b4)
    OUT <- rbind(OUT, save)
  }
  plot(0, 0, xlim = c(0.6, 4.4), ylim = c(0, 1), xaxt = 'n', main = demo, ylab = '', xlab = '', yaxt = 'n', bty = 'n')
  axis(1, at = c(1:4), lwd = lwd, cex.axis = cex.axis)
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c('0.0','','0.5','','1.0'), lwd = lwd, cex.axis = cex.axis)
  x <- 1
  for(c in 2:5){
    points(jitter(rep(x, 50), amount = .15), OUT[,c], pch = 16, col = alpha('springgreen4', pt.alph), cex = pt.cex)
    lines(c(x - line.len, x + line.len), c(median(OUT[,c]), median(OUT[,c])), lwd = 3, col = 'springgreen4')
    # lines(c(x - line.len, x + line.len), c(median(OUT[,c]), median(OUT[,c])), lwd = 3, col = 'springgreen3')
    x <- x+1
  }
}
dev.off()

