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
          save <- c(samp, s, e, (e-s+1))
          OUT1 <- rbind(OUT1, save)
        }
        if(h == nrow(het.sites) & nrow(het.sites)!= 1){
          s <- het.sites[h-1, 1] + 1
          e <- het.sites[h, 1] - 1
          save <- c(samp, s, e, (e-s+1))
          OUT1 <- rbind(OUT1, save)
          s <- het.sites[h, 1] + 1
          e <- c.len
          save <- c(samp, s, e, (e-s+1))
          OUT1 <- rbind(OUT1, save)
        }
        if(h != 1 & h != nrow(het.sites) & nrow(het.sites)!= 1){
          s <- het.sites[h-1, 1] + 1
          e <- het.sites[h, 1] - 1
          save <- c(samp, s, e, (e-s+1))
          OUT1 <- rbind(OUT1, save)
        }
        if(nrow(het.sites) == 1){
          s <- 1
          e <- het.sites[h, 1] - 1
          save <- c(samp, s, e, (e-s+1))
          OUT1 <- rbind(OUT1, save)
          
          s <- het.sites[h, 1] + 1
          e <- c.len
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

### Format histograms for manuscript figures
pdf('/Users/Avril/Desktop/scenario_froh.pdf', width = 5, height = 5)
bottle <- read.table('bottle_true_roh_coords.txt', header = FALSE)
colnames(bottle) <- c('id','start','end','length')
bottle <- bottle[bottle$length >= 100e3,]
SAVE <- NULL
for(i in unique(bottle$id)){
  froh <- sum(bottle[bottle$id == i, 'length'])/c.len
  SAVE <- c(SAVE, froh)
}
hist(SAVE, xlim = c(0,1), xlab = '', main = 'bottle', ylim = c(0,25), breaks = 10) ## 19

decline <- read.table('decline_true_roh_coords.txt', header = FALSE)
colnames(decline) <- c('id','start','end','length')
decline <- decline[decline$length >= 100e3,]
SAVE <- NULL
for(i in unique(decline$id)){
  froh <- sum(decline[decline$id == i, 'length'])/c.len
  SAVE <- c(SAVE, froh)
}
hist(SAVE, xlim = c(0,1), xlab = '', main = 'decline', ylim = c(0,25), breaks = 11) ## 19

small <- read.table('small_true_roh_coords.txt', header = FALSE)
colnames(small) <- c('id','start','end','length')
small <- small[small$length >= 100e3,]
SAVE <- NULL
for(i in unique(small$id)){
  froh <- sum(small[small$id == i, 'length'])/c.len
  SAVE <- c(SAVE, froh)
}
hist(SAVE, xlim = c(0,1), xlab = '', main = 'small', ylim = c(0,25), breaks = 8) ## 9

large.1000 <- read.table('large-1000_true_roh_coords.txt', header = FALSE)
colnames(large.1000) <- c('id','start','end','length')
large.1000 <- large.1000[large.1000$length >= 100e3,]
SAVE <- NULL
for(i in unique(large.1000$id)){
  froh <- sum(large.1000[large.1000$id == i, 'length'])/c.len
  SAVE <- c(SAVE, froh)
}
hist(SAVE, xlim = c(0,1), xlab = '', main = 'large.1000', ylim = c(0,25), breaks = 7) ## 8
dev.off()

b.1 <- 100e3
b.2 <- 500e3
b.3 <- 1e6
b.4 <- 2e6

pdf('/Users/Avril/Desktop/scenario_froh_length_bins.pdf', width = 4, height = 5)
pt.cex = 0.8

bottle <- read.table('bottle_true_roh_coords.txt', header = FALSE)
colnames(bottle) <- c('id','start','end','length')
bottle <- bottle[bottle$length >= 100e3,]
OUT <- NULL
for(i in unique(bottle$id)){
  sub <- bottle[bottle$id == i,]
  b1 <- sum(sub[sub$length >= b.1 & sub$length < b.2, 'length'])/c.len
  b2 <- sum(sub[sub$length >= b.2 & sub$length < b.3, 'length'])/c.len
  b3 <- sum(sub[sub$length >= b.3 & sub$length < b.4, 'length'])/c.len
  b4 <- sum(sub[sub$length >= b.4, 'length'])/c.len
  
  save <- c(i, b1, b2, b3, b4)
  OUT <- rbind(OUT, save)
}
plot(0, 0, xlim = c(0.75, 4.25), ylim = c(0, 1), xaxt = 'n', main = 'bottle', ylab = '', xlab = '', yaxt = 'n')
  axis(1, at = c(1:4))
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1))
  points(jitter(rep(1, 50), amount = .15), OUT[,2], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 5, col = 'springgreen4')
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(2, 50), amount = .15), OUT[,3], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 5, col = 'springgreen4')
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(3, 50), amount = .15), OUT[,4], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 5, col = 'springgreen4')
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(4, 50), amount = .15), OUT[,5], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 5, col = 'springgreen4')
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 3, col = 'springgreen3')

decline <- read.table('decline_true_roh_coords.txt', header = FALSE)
colnames(decline) <- c('id','start','end','length')
decline <- decline[decline$length >= 100e3,]
OUT <- NULL
for(i in unique(decline$id)){
  sub <- decline[decline$id == i,]
  b1 <- sum(sub[sub$length >= b.1 & sub$length < b.2, 'length'])/c.len
  b2 <- sum(sub[sub$length >= b.2 & sub$length < b.3, 'length'])/c.len
  b3 <- sum(sub[sub$length >= b.3 & sub$length < b.4, 'length'])/c.len
  b4 <- sum(sub[sub$length >= b.4, 'length'])/c.len
  
  save <- c(i, b1, b2, b3, b4)
  OUT <- rbind(OUT, save)
}
plot(0, 0, xlim = c(0.75, 4.25), ylim = c(0, 1), xaxt = 'n', main = 'decline', ylab = '', xlab = '', yaxt = 'n')
  axis(1, at = c(1:4))
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1))
  points(jitter(rep(1, 50), amount = .15), OUT[,2], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 5, col = 'springgreen4')
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(2, 50), amount = .15), OUT[,3], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 5, col = 'springgreen4')
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(3, 50), amount = .15), OUT[,4], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 5, col = 'springgreen4')
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(4, 50), amount = .15), OUT[,5], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 5, col = 'springgreen4')
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 3, col = 'springgreen3')

small <- read.table('small_true_roh_coords.txt', header = FALSE)
colnames(small) <- c('id','start','end','length')
small <- small[small$length >= 100e3,]
OUT <- NULL
for(i in unique(small$id)){
  sub <- small[small$id == i,]
  b1 <- sum(sub[sub$length >= b.1 & sub$length < b.2, 'length'])/c.len
  b2 <- sum(sub[sub$length >= b.2 & sub$length < b.3, 'length'])/c.len
  b3 <- sum(sub[sub$length >= b.3 & sub$length < b.4, 'length'])/c.len
  b4 <- sum(sub[sub$length >= b.4, 'length'])/c.len
  
  save <- c(i, b1, b2, b3, b4)
  OUT <- rbind(OUT, save)
}
plot(0, 0, xlim = c(0.75, 4.25), ylim = c(0, 1), xaxt = 'n', main = 'small', ylab = '', xlab = '', yaxt = 'n')
  axis(1, at = c(1:4))
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1))
  points(jitter(rep(1, 50), amount = .15), OUT[,2], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 5, col = 'springgreen4')
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(2, 50), amount = .15), OUT[,3], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 5, col = 'springgreen4')
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(3, 50), amount = .15), OUT[,4], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 5, col = 'springgreen4')
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(4, 50), amount = .15), OUT[,5], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 5, col = 'springgreen4')
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 3, col = 'springgreen3')

large.1000 <- read.table('large-1000_true_roh_coords.txt', header = FALSE)
colnames(large.1000) <- c('id','start','end','length')
large.1000 <- large.1000[large.1000$length >= 100e3,]
OUT <- NULL
for(i in unique(large.1000$id)){
  sub <- large.1000[large.1000$id == i,]
  b1 <- sum(sub[sub$length >= b.1 & sub$length < b.2, 'length'])/c.len
  b2 <- sum(sub[sub$length >= b.2 & sub$length < b.3, 'length'])/c.len
  b3 <- sum(sub[sub$length >= b.3 & sub$length < b.4, 'length'])/c.len
  b4 <- sum(sub[sub$length >= b.4, 'length'])/c.len
  
  save <- c(i, b1, b2, b3, b4)
  OUT <- rbind(OUT, save)
}
plot(0, 0, xlim = c(0.75, 4.25), ylim = c(0, 1), xaxt = 'n', main = 'large.1000', ylab = '', xlab = '', yaxt = 'n')
  axis(1, at = c(1:4))
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1))
  points(jitter(rep(1, 50), amount = .15), OUT[,2], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 5, col = 'springgreen4')
  lines(c(0.7, 1.3), c(median(OUT[,2]), median(OUT[,2])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(2, 50), amount = .15), OUT[,3], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 5, col = 'springgreen4')
  lines(c(1.7, 2.3), c(median(OUT[,3]), median(OUT[,3])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(3, 50), amount = .15), OUT[,4], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 5, col = 'springgreen4')
  lines(c(2.7, 3.3), c(median(OUT[,4]), median(OUT[,4])), lwd = 3, col = 'springgreen3')
  points(jitter(rep(4, 50), amount = .15), OUT[,5], pch = 19, col = alpha('springgreen4', 0.6), cex = pt.cex)
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 5, col = 'springgreen4')
  lines(c(3.7, 4.3), c(median(OUT[,5]), median(OUT[,5])), lwd = 3, col = 'springgreen3')

dev.off()
