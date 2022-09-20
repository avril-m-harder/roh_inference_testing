##### Code for analyzing empirical data
library(scales)
library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggridges)
theme_set(theme_bw(16))
`%notin%` <- Negate(`%in%`)

### set colors
library(ghibli)
gt.col <- ghibli_palette('PonyoMedium')[4]
lt.gt.col <- ghibli_palette('PonyoLight')[4]
pl.col <- ghibli_palette('PonyoMedium')[2]
lt.pl.col <- ghibli_palette('PonyoLight')[2]
plink.col <- ghibli_palette('PonyoMedium')[3]
lt.plink.col <- ghibli_palette('PonyoLight')[3]
min.covg.col <- ghibli_palette('PonyoMedium')[6]
max.covg.col <- ghibli_palette('PonyoMedium')[1]

setwd('/Users/Avril/Desktop/roh_param_project/empirical_plotting_summarizing/')

##### 1A. Read in data and summary files #####
## tas dev chrom data
chroms <- read.table('tasdev_assembly_data/mSarHar1.11_autosomes.txt', sep = '\t', header = TRUE)
chroms <- chroms[,c(5,9)]
colnames(chroms) <- c('chrom','length')
chroms$cum.len <- cumsum(as.numeric(chroms$length))
tot.len <- sum(chroms$length)

## Sample and chrom keys
samp.key <- read.csv('../empirical_qc/sample_id_key.csv')
chrom.key <- read.csv('../empirical_qc/chrom_key.csv')

## BCFtools/ROH results
gt.res <- read.csv('bcftools_data/tasdev_GT_allcovg_levels.csv')
gt.res <- merge(gt.res, chrom.key, by.x = 'chrom', by.y = 'chrom.name')
pl.res <- read.csv('bcftools_data/tasdev_PL_allcovg_levels.csv')
pl.res <- merge(pl.res, chrom.key, by.x = 'chrom', by.y = 'chrom.name')

## settled on default settings
plink.res <- read.csv('plink_results_round1/PLINK_all_coordinates_DEFAULT_SETTINGS.csv')
plink.res$length <- plink.res$end - plink.res$start + 1

## Individual f(ROH) values
froh.results <- read.csv('combined_empirical_data/individual_froh_stats.csv')


##### 1B. Plot histograms of ROH length distributions #####
alph <- 0.4
pdf('../manuscript/r_scripts_AMH/figures_output/empirical/ROH_length_histograms.pdf', width = 6, height = 5)

## called BCFtools Genotypes
hist(gt.res$length, breaks = 50, xlab = 'Called ROH length (Mb)', main = 'BCFtools Genotypes - empirical data', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
hist(gt.res$length, breaks = 100, xlab = 'Called ROH length', main = 'BCFtools Genotypes - empirical data', xlim = c(1e5, 1e6))

## called BCFtools Likelihoods
hist(pl.res$length, breaks = 25, xlab = 'Called ROH length (Mb)', main = 'BCFtools Likelihoods - empirical data', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6), labels = c(0, 1, 2, 3, 4, 5, 6))
hist(pl.res$length, breaks = 100, xlab = 'Called ROH length', main = 'BCFtools Likelihoods - empirical data', xlim = c(1e5, 1e6))

## called PLINK
hist(plink.res$length, breaks = 50, xlab = 'Called ROH length (Mb)', main = 'PLINK - empirical data', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7))
hist(plink.res$length, breaks = 100, xlab = 'Called ROH length', main = 'PLINK - empirical data', xlim = c(1e5, 1e6))

## combined
hist(gt.res$length, breaks = 100, xlab = 'True ROH length (kb)', main = 'True empirical data', xlim = c(1e5, 1e6), 
     border = 'transparent', col = 'transparent', xaxt = 'n')
  axis(1, at = c(0, 2.5e5, 5e5, 7.5e5, 1e6), labels = c(0, 250, 500, 750, 1000))
  hist(plink.res$length, breaks = 100, add = TRUE, col = alpha(plink.col, alph))
  hist(pl.res$length, breaks = 100, add = TRUE, col = alpha(pl.col, alph))
  hist(gt.res$length, breaks = 200, add = TRUE, col = alpha(gt.col, alph))
  
plot(density(pl.res$length), xlim = c(1e5, 7.5e6), xaxt = 'n', xlab = 'ROH length (Mb)', col = 'transparent')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
  polygon(density(pl.res$length), col = alpha(pl.col, alph), border = pl.col)
  polygon(density(plink.res$length), col = alpha(plink.col, alph), border = plink.col)
  polygon(density(gt.res$length), col = alpha(gt.col, alph), border = gt.col)
  legend('topright', legend = c('Genotypes', 'Likelihoods','PLINK'), border = c(gt.col, pl.col, plink.col),
         fill = c(alpha(gt.col, alph), alpha(pl.col, alph), alpha(plink.col, alph)), bty = 'n')
  
plot(density(pl.res$length), xlim = c(1e5, 1e6), xaxt = 'n', xlab = 'ROH length (kb)', col = 'transparent')
  axis(1, at = c(0, 2.5e5, 5e5, 7.5e5, 1e6), labels = c(0, 250, 500, 750, 1000))
  polygon(density(pl.res$length), col = alpha(pl.col, alph), border = pl.col)
  polygon(density(plink.res$length), col = alpha(plink.col, alph), border = plink.col)
  polygon(density(gt.res$length), col = alpha(gt.col, alph), border = gt.col)
  legend('topright', legend = c('Genotypes', 'Likelihoods','PLINK'), border = c(gt.col, pl.col, plink.col),
         fill = c(alpha(gt.col, alph), alpha(pl.col, alph), alpha(plink.col, alph)), bty = 'n')

plot(density(pl.res$length), xlim = c(1e5, 5e5), xaxt = 'n', xlab = 'ROH length (kb)', col = 'transparent')
  axis(1, at = c(0, 1e5, 2e5, 3e5, 4e5, 5e5), labels = c(0, 100, 200, 300, 400, 500))
  polygon(density(pl.res$length), col = alpha(pl.col, alph), border = pl.col) 
  polygon(density(plink.res$length), col = alpha(plink.col, alph), border = plink.col)
  polygon(density(gt.res$length), col = alpha(gt.col, alph), border = gt.col)
  legend('topright', legend = c('Genotypes', 'Likelihoods','PLINK'), border = c(gt.col, pl.col, plink.col),
         fill = c(alpha(gt.col, alph), alpha(pl.col, alph), alpha(plink.col, alph)), bty = 'n')

dev.off()

##### 1C. Plot cumulative f(ROH) for each sample #####
# pdf('../manuscript/r_scripts_AMH/figures_output/empirical/cumulative_fROH_figures_empirical.pdf', width = 8, height = 6)

pdf('../manuscript/r_scripts_AMH/figures_output/empirical/cumulative_fROH_figures_empirical.pdf', width = 20, height = 12)
par(mfrow = c(3, 5))

## Genotypes
for(c in unique(sort(gt.res$covg))){
  plot(0, 0, col = 'transparent', xlim = c(1e5, 7.5e6), ylim = c(0, max(froh.results[froh.results$covg == c, 'gt.froh'])), 
       main = paste0('Empirical: BCFtools Genotypes - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
  for(i in unique(gt.res$id)){
    sub <- gt.res[gt.res$id == i & gt.res$covg == c,]
    sub <- sub[order(sub$length),]
    sub$cumlen <- cumsum(sub$length)
    sub$cumfroh <- sub$cumlen / tot.len
    lines(sub$length, sub$cumfroh, col = alpha(gt.col, 0.6))
  }
}

## Likelihoods
for(c in unique(sort(pl.res$covg))){
  plot(0, 0, col = 'transparent', xlim = c(1e5, 7.5e6), ylim = c(0, max(froh.results[froh.results$covg == c, 'pl.froh'])), 
       main = paste0('Empirical: BCFtools Likelihoods - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
  for(i in unique(pl.res$id)){
    sub <- pl.res[pl.res$id == i & pl.res$covg == c,]
    sub <- sub[order(sub$length),]
    sub$cumlen <- cumsum(sub$length)
    sub$cumfroh <- sub$cumlen / tot.len
    lines(sub$length, sub$cumfroh, col = alpha(pl.col, 0.6))
  }
}

## PLINK
for(c in unique(sort(plink.res$covg))){
  plot(0, 0, col = 'transparent', xlim = c(1e5, 7.5e6), ylim = c(0, max(froh.results[froh.results$covg == c, 'plink.froh'])), 
       main = paste0('Empirical: PLINK - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
  for(i in unique(plink.res$id)){
    sub <- plink.res[plink.res$id == i & plink.res$covg == c,]
    sub <- sub[order(sub$length),]
    sub$cumlen <- cumsum(sub$length)
    sub$cumfroh <- sub$cumlen / tot.len
    lines(sub$length, sub$cumfroh, col = alpha(plink.col, 0.6))
  }
}

dev.off()

# pdf('../manuscript/r_scripts_AMH/figures_output/empirical/cumulative_fROH_figures_empirical_uniformyaxis.pdf', width = 8, height = 6)

pdf('../manuscript/r_scripts_AMH/figures_output/empirical/cumulative_fROH_figures_empirical_uniformyaxis.pdf', width = 20, height = 12)
par(mfrow = c(3,5))

pt.cex <- 0.2
alph <- 0.4

## Genotypes
for(c in unique(sort(gt.res$covg))){
  plot(0, 0, col = 'transparent', xlim = c(1e5, 7.5e6), ylim = c(0, 0.42), 
       main = paste0('Empirical: BCFtools Genotypes - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
  for(i in unique(gt.res$id)){
    sub <- gt.res[gt.res$id == i & gt.res$covg == c,]
    sub <- sub[order(sub$length),]
    sub$cumlen <- cumsum(sub$length)
    sub$cumfroh <- sub$cumlen / tot.len
    lines(sub$length, sub$cumfroh, col = alpha(gt.col, alph))
    points(sub$length, sub$cumfroh, col = gt.col, cex = pt.cex)
  }
}

## Likelihoods
for(c in unique(sort(pl.res$covg))){
  plot(0, 0, col = 'transparent', xlim = c(1e5, 7.5e6), ylim = c(0, 0.42), 
       main = paste0('Empirical: BCFtools Likelihoods - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
  for(i in unique(pl.res$id)){
    sub <- pl.res[pl.res$id == i & pl.res$covg == c,]
    sub <- sub[order(sub$length),]
    sub$cumlen <- cumsum(sub$length)
    sub$cumfroh <- sub$cumlen / tot.len
    lines(sub$length, sub$cumfroh, col = alpha(pl.col, alph))
    points(sub$length, sub$cumfroh, col = pl.col, cex = pt.cex)
  }
}

## PLINK
for(c in unique(sort(plink.res$covg))){
  plot(0, 0, col = 'transparent', xlim = c(1e5, 7.5e6), ylim = c(0, 0.42), 
       main = paste0('Empirical: PLINK - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
  axis(1, at = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6), labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
  for(i in unique(plink.res$id)){
    sub <- plink.res[plink.res$id == i & plink.res$covg == c,]
    sub <- sub[order(sub$length),]
    sub$cumlen <- cumsum(sub$length)
    sub$cumfroh <- sub$cumlen / tot.len
    lines(sub$length, sub$cumfroh, col = alpha(plink.col, alph))
    points(sub$length, sub$cumfroh, col = plink.col, cex = pt.cex)
  }
}

dev.off()


##### 2. Plotting f(ROH) by length bins, individual lines #####
b1 <- 5e5
b2 <- 1e6
b3 <- 2e6

ymin <- 0
ymax <- 0.15
xmin <- 0.85
xmax <- 5.15
alph <- 0.3

k <- 1
while(k == 1){
  pdf(paste0('../manuscript/r_scripts_AMH/figures_output/empirical/fROH_by_length_bins_indivlines.pdf'), width = 12, height = 12)
  par(mfrow = c(3,4))
  
  ## GT
  OUT <- NULL
  for(i in unique(gt.res$id)){
    for(c in unique(gt.res$covg)){
      temp <- gt.res[gt.res$id == i & gt.res$covg == c,]
      save <- c(i, c, sum(temp$length)/tot.len, nrow(temp), 
                sum(temp[temp$length < b1, 'length'])/tot.len, nrow(temp[temp$length < b1,]),
                sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/tot.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/tot.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                sum(temp[temp$length >= b3, 'length'])/tot.len, nrow(temp[temp$length >= b3,]))
      OUT <- rbind(OUT, save)
    }
  }
  gt.frohs <- as.data.frame(OUT)
  colnames(gt.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
  for(c in 2:ncol(gt.frohs)){
    gt.frohs[,c] <- as.numeric(gt.frohs[,c])
  }
  
  ## GT plot
  # par(mfrow = c(2,4))
  ## bin1
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Genotypes only: Short ROHs', xlab = 'Coverage', ylab = 'f(ROH)',
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin1.froh, pch = 19, col = gt.col)
    lines(c(1:5), temp$bin1.froh, col = alpha(gt.col, alph))
  }
  ## bin 2
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax),  
       xaxt = 'n', main = 'Genotypes only: Intermediate ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin2.froh, pch = 19, col = gt.col)
    lines(c(1:5), temp$bin2.froh, col = alpha(gt.col, alph))
  }
  ## bin 3
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Genotypes only: Long ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin3.froh, pch = 19, col = gt.col)
    lines(c(1:5), temp$bin3.froh, col = alpha(gt.col, alph))
  }
  ## bin 4
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Genotypes only: Longer ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin4.froh, pch = 19, col = gt.col)
    lines(c(1:5), temp$bin4.froh, col = alpha(gt.col, alph))
  }
  
  ## PL
  OUT <- NULL
  for(i in unique(pl.res$id)){
    for(c in unique(pl.res$covg)){
      temp <- pl.res[pl.res$id == i & pl.res$covg == c,]
      save <- c(i, c, sum(temp$length)/tot.len, nrow(temp), 
                sum(temp[temp$length < b1, 'length'])/tot.len, nrow(temp[temp$length < b1,]),
                sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/tot.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/tot.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                sum(temp[temp$length >= b3, 'length'])/tot.len, nrow(temp[temp$length >= b3,]))
      OUT <- rbind(OUT, save)
    }
  }
  pl.frohs <- as.data.frame(OUT)
  colnames(pl.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
  for(c in 2:ncol(pl.frohs)){
    pl.frohs[,c] <- as.numeric(pl.frohs[,c])
  }
  
  ## PL plot
  ## bin1
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Genotype likelihoods: Short ROHs', xlab = 'Coverage', ylab = 'f(ROH)',
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(pl.frohs$id)){
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin1.froh, pch = 19, col = pl.col)
    lines(c(1:5), temp$bin1.froh, col = alpha(pl.col, alph))
  }
  ## bin 2
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax),  
       xaxt = 'n', main = 'Genotype likelihoods: Intermediate ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(pl.frohs$id)){
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin2.froh, pch = 19, col = pl.col)
    lines(c(1:5), temp$bin2.froh, col = alpha(pl.col, alph))
  }
  ## bin 3
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Genotype likelihoods: Long ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(pl.frohs$id)){
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin3.froh, pch = 19, col = pl.col)
    lines(c(1:5), temp$bin3.froh, col = alpha(pl.col, alph))
  }
  ## bin 4
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Genotype likelihoods: Longer ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(pl.frohs$id)){
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin4.froh, pch = 19, col = pl.col)
    lines(c(1:5), temp$bin4.froh, col = alpha(pl.col, alph))
  }
  
  ## PLINK
  OUT <- NULL
  for(i in unique(plink.res$id)){
    for(c in unique(plink.res$covg)){
      temp <- plink.res[plink.res$id == i & plink.res$covg == c,]
      save <- c(i, c, sum(temp$length)/tot.len, nrow(temp), 
                sum(temp[temp$length < b1, 'length'])/tot.len, nrow(temp[temp$length < b1,]),
                sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/tot.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/tot.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                sum(temp[temp$length >= b3, 'length'])/tot.len, nrow(temp[temp$length >= b3,]))
      OUT <- rbind(OUT, save)
    }
  }
  plink.frohs <- as.data.frame(OUT)
  colnames(plink.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
  for(c in 2:ncol(plink.frohs)){
    plink.frohs[,c] <- as.numeric(plink.frohs[,c])
  }
  
  ## PLINK plot
  ## bin1
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'PLINK: Short ROHs', xlab = 'Coverage', ylab = 'f(ROH)',
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(plink.frohs$id)){
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin1.froh, pch = 19, col = plink.col)
    lines(c(1:5), temp$bin1.froh, col = alpha(plink.col, alph))
  }
  ## bin 2
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax),  
       xaxt = 'n', main = 'PLINK: Intermediate ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(plink.frohs$id)){
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin2.froh, pch = 19, col = plink.col)
    lines(c(1:5), temp$bin2.froh, col = alpha(plink.col, alph))
  }
  ## bin 3
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'PLINK: Long ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(plink.frohs$id)){
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin3.froh, pch = 19, col = plink.col)
    lines(c(1:5), temp$bin3.froh, col = alpha(plink.col, alph))
  }
  ## bin 4
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'PLINK: Longer ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(plink.frohs$id)){
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    points(c(1:5), temp$bin4.froh, pch = 19, col = plink.col)
    lines(c(1:5), temp$bin4.froh, col = alpha(plink.col, alph))
  }
  
  k <- k+1
  dev.off()
}

## Calculate mean difference in f(ROH) between 5X and 10X for all ROH length bins, combining methods
OUT <- NULL
for(i in unique(plink.frohs$id)){
  sub.plink <- plink.frohs[plink.frohs$id == i,]
  sub.gt <- gt.frohs[gt.frohs$id == i,]
  sub.pl <- pl.frohs[pl.frohs$id == i,]
  save <- c(sub.pl[sub.pl$covg == 5, 'bin1.froh'] - sub.pl[sub.pl$covg == 10, 'bin1.froh'],
            sub.pl[sub.pl$covg == 5, 'bin2.froh'] - sub.pl[sub.pl$covg == 10, 'bin2.froh'],
            sub.pl[sub.pl$covg == 5, 'bin3.froh'] - sub.pl[sub.pl$covg == 10, 'bin3.froh'],
            sub.pl[sub.pl$covg == 5, 'bin4.froh'] - sub.pl[sub.pl$covg == 10, 'bin4.froh'],
            sub.gt[sub.gt$covg == 5, 'bin1.froh'] - sub.gt[sub.gt$covg == 10, 'bin1.froh'],
            sub.gt[sub.gt$covg == 5, 'bin2.froh'] - sub.gt[sub.gt$covg == 10, 'bin2.froh'],
            sub.gt[sub.gt$covg == 5, 'bin3.froh'] - sub.gt[sub.gt$covg == 10, 'bin3.froh'],
            sub.gt[sub.gt$covg == 5, 'bin4.froh'] - sub.gt[sub.gt$covg == 10, 'bin4.froh'],
            sub.plink[sub.plink$covg == 5, 'bin1.froh'] - sub.plink[sub.plink$covg == 10, 'bin1.froh'],
            sub.plink[sub.plink$covg == 5, 'bin2.froh'] - sub.plink[sub.plink$covg == 10, 'bin2.froh'],
            sub.plink[sub.plink$covg == 5, 'bin3.froh'] - sub.plink[sub.plink$covg == 10, 'bin3.froh'],
            sub.plink[sub.plink$covg == 5, 'bin4.froh'] - sub.plink[sub.plink$covg == 10, 'bin4.froh'])
  OUT <- rbind(OUT, save)
}
cats <- c('PL short','PL intermediate','PL long','PL very long',
          'GT short','GT intermediate','GT long','GT very long',
          'PLINK short','PLINK intermediate','PLINK long','PLINK very long')
colnames(OUT) <- cats
mean(colMeans(OUT)[c(1,5,9)]) ## short
mean(colMeans(OUT)[c(2,6,10)]) ## intermediate
mean(colMeans(OUT)[c(3,7,11)]) ## long
mean(colMeans(OUT)[c(4,8,12)]) ## very long

##### >>> 2A. Plotting f(ROH) by length bins, 1 plot per bin #####
b1 <- 5e5
b2 <- 1e6
b3 <- 2e6

ymin <- 0
ymax <- 0.15
xmin <- 0.85
xmax <- 5.15
bg.alph <- 0.25
lwd <- 2
offset <- 0.15
err.width <- 0.05
txt.size <- 1.25

k <- 1
while(k == 1){
  pdf(paste0('../manuscript/r_scripts_AMH/figures_output/empirical/fROH_by_length_bins_indivlines.pdf'), width = 5, height = 6)
  par(mfrow = c(1,1), mar = c(5.1, 4.6, 4.1, 2.1))
  
  ## GT
  OUT <- NULL
  for(i in unique(gt.res$id)){
    for(c in unique(gt.res$covg)){
      temp <- gt.res[gt.res$id == i & gt.res$covg == c,]
      save <- c(i, c, sum(temp$length)/tot.len, nrow(temp), 
                sum(temp[temp$length < b1, 'length'])/tot.len, nrow(temp[temp$length < b1,]),
                sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/tot.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/tot.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                sum(temp[temp$length >= b3, 'length'])/tot.len, nrow(temp[temp$length >= b3,]))
      OUT <- rbind(OUT, save)
    }
  }
  gt.frohs <- as.data.frame(OUT)
  colnames(gt.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
  for(c in 2:ncol(gt.frohs)){
    gt.frohs[,c] <- as.numeric(gt.frohs[,c])
  }
  
  ## PL
  OUT <- NULL
  for(i in unique(pl.res$id)){
    for(c in unique(pl.res$covg)){
      temp <- pl.res[pl.res$id == i & pl.res$covg == c,]
      save <- c(i, c, sum(temp$length)/tot.len, nrow(temp), 
                sum(temp[temp$length < b1, 'length'])/tot.len, nrow(temp[temp$length < b1,]),
                sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/tot.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/tot.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                sum(temp[temp$length >= b3, 'length'])/tot.len, nrow(temp[temp$length >= b3,]))
      OUT <- rbind(OUT, save)
    }
  }
  pl.frohs <- as.data.frame(OUT)
  colnames(pl.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
  for(c in 2:ncol(pl.frohs)){
    pl.frohs[,c] <- as.numeric(pl.frohs[,c])
  }
  
  ## PLINK
  OUT <- NULL
  for(i in unique(plink.res$id)){
    for(c in unique(plink.res$covg)){
      temp <- plink.res[plink.res$id == i & plink.res$covg == c,]
      save <- c(i, c, sum(temp$length)/tot.len, nrow(temp), 
                sum(temp[temp$length < b1, 'length'])/tot.len, nrow(temp[temp$length < b1,]),
                sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/tot.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/tot.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                sum(temp[temp$length >= b3, 'length'])/tot.len, nrow(temp[temp$length >= b3,]))
      OUT <- rbind(OUT, save)
    }
  }
  plink.frohs <- as.data.frame(OUT)
  colnames(plink.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
  for(c in 2:ncol(plink.frohs)){
    plink.frohs[,c] <- as.numeric(plink.frohs[,c])
  }
  
  ## Bin 1 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Short ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin1.froh, col = alpha(gt.col, bg.alph))
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin1.froh, col = alpha(pl.col, bg.alph))
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin1.froh, col = alpha(plink.col, bg.alph))
  }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin1.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin1.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin1.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin1.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin1.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin1.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')
  write.csv(OUT, '../manuscript/r_scripts_AMH/tables_output/empirical_bin1.csv')
  
  ## SEs
  # lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  # points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  # arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se),
  #        y1 = (OUT$gt.mean + OUT$gt.se),
  #        lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  # points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  # arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se),
  #        y1 = (OUT$pl.mean + OUT$pl.se),
  #        lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  # points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  # arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se),
  #        y1 = (OUT$plink.mean + OUT$plink.se),
  #        lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  legend('bottomleft', pch = 19, lwd = lwd, legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'), bty = 'n', inset = 0.02, 
         col = c(gt.col, pl.col, plink.col), cex = txt.size)
  
  ## Bin 2 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Intermediate ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin2.froh, col = alpha(gt.col, bg.alph))
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin2.froh, col = alpha(pl.col, bg.alph))
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin2.froh, col = alpha(plink.col, bg.alph))
  }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin2.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin2.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin2.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin2.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin2.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin2.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')
  write.csv(OUT, '../manuscript/r_scripts_AMH/tables_output/empirical_bin2.csv')
  
  ## SEs
  # lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  # points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  # arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se),
  #        y1 = (OUT$gt.mean + OUT$gt.se),
  #        lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  # points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  # arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se),
  #        y1 = (OUT$pl.mean + OUT$pl.se),
  #        lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  # points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  # arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se),
  #        y1 = (OUT$plink.mean + OUT$plink.se),
  #        lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## Bin 3 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Long ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin3.froh, col = alpha(gt.col, bg.alph))
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin3.froh, col = alpha(pl.col, bg.alph))
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin3.froh, col = alpha(plink.col, bg.alph))
  }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin3.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin3.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin3.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin3.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin3.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin3.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')
  write.csv(OUT, '../manuscript/r_scripts_AMH/tables_output/empirical_bin3.csv')
  
  ## SEs  
  # lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  # points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  # arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se),
  #        y1 = (OUT$gt.mean + OUT$gt.se),
  #        lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  # points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  # arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se),
  #        y1 = (OUT$pl.mean + OUT$pl.se),
  #        lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  # points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  # arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se),
  #        y1 = (OUT$plink.mean + OUT$plink.se),
  #        lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## Bin 4 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Very long ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  for(i in unique(gt.frohs$id)){
    temp <- gt.frohs[gt.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin4.froh, col = alpha(gt.col, bg.alph))
    temp <- pl.frohs[pl.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin4.froh, col = alpha(pl.col, bg.alph))
    temp <- plink.frohs[plink.frohs$id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), temp$bin4.froh, col = alpha(plink.col, bg.alph))
  }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin4.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin4.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin4.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin4.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin4.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin4.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')
  write.csv(OUT, '../manuscript/r_scripts_AMH/tables_output/empirical_bin4.csv')
  
  ## SEs
  # lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  # points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  # arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se),
  #        y1 = (OUT$gt.mean + OUT$gt.se),
  #        lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  # points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  # arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se),
  #        y1 = (OUT$pl.mean + OUT$pl.se),
  #        lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  # 
  # lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  # points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  # arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se),
  #        y1 = (OUT$plink.mean + OUT$plink.se),
  #        lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)

  # legend('topright', pch = 19, lwd = lwd, legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'), bty = 'n', inset = c(0, 0.07), 
  #        col = c(gt.col, pl.col, plink.col), cex = txt.size)

  dev.off()
  
  ### plots for stats comps of methods in SI
  pdf(paste0('../manuscript/r_scripts_AMH/figures_output/empirical/fROH_by_length_bins_SI_method_comp_plot.pdf'), width = 5, height = 6)
  par(mfrow = c(1,1), mar = c(5.1, 4.6, 4.1, 2.1))
  
  offset <- 0.2
  
  ## Bin 1 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Short ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  # for(i in unique(gt.frohs$id)){
  #   temp <- gt.frohs[gt.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin1.froh, col = alpha(gt.col, bg.alph))
  #   temp <- pl.frohs[pl.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin1.froh, col = alpha(pl.col, bg.alph))
  #   temp <- plink.frohs[plink.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin1.froh, col = alpha(plink.col, bg.alph))
  # }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin1.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin1.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin1.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin1.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin1.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin1.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')

  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  legend('bottomleft', pch = 19, lwd = lwd, legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'), bty = 'n', inset = 0.02, 
         col = c(gt.col, pl.col, plink.col), cex = txt.size)
  
  ## Bin 2 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Intermediate ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  # for(i in unique(gt.frohs$id)){
  #   temp <- gt.frohs[gt.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin2.froh, col = alpha(gt.col, bg.alph))
  #   temp <- pl.frohs[pl.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin2.froh, col = alpha(pl.col, bg.alph))
  #   temp <- plink.frohs[plink.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin2.froh, col = alpha(plink.col, bg.alph))
  # }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin2.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin2.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin2.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin2.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin2.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin2.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')
  
  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## Bin 3 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Long ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  # for(i in unique(gt.frohs$id)){
  #   temp <- gt.frohs[gt.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin3.froh, col = alpha(gt.col, bg.alph))
  #   temp <- pl.frohs[pl.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin3.froh, col = alpha(pl.col, bg.alph))
  #   temp <- plink.frohs[plink.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin3.froh, col = alpha(plink.col, bg.alph))
  # }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin3.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin3.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin3.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin3.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin3.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin3.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')
  
  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  ## Bin 4 plot
  plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
       xaxt = 'n', main = 'Very long ROHs', xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
       cex.axis = 1.25, cex.lab = 1.25)
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = 1.25)
  x <- 1
  # for(i in unique(gt.frohs$id)){
  #   temp <- gt.frohs[gt.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin4.froh, col = alpha(gt.col, bg.alph))
  #   temp <- pl.frohs[pl.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin4.froh, col = alpha(pl.col, bg.alph))
  #   temp <- plink.frohs[plink.frohs$id == i,]
  #   temp <- temp[order(temp$covg),]
  #   lines(c(1:5), temp$bin4.froh, col = alpha(plink.col, bg.alph))
  # }
  OUT <- NULL
  for(c in sort(unique(gt.frohs$covg))){
    save <- c(c, mean(gt.frohs[gt.frohs$covg == c, 'bin4.froh']), sd(gt.frohs[gt.frohs$covg == c, 'bin4.froh'])/sqrt(15),
              mean(pl.frohs[pl.frohs$covg == c, 'bin4.froh']), sd(pl.frohs[pl.frohs$covg == c, 'bin4.froh'])/sqrt(15),
              mean(plink.frohs[plink.frohs$covg == c, 'bin4.froh']), sd(plink.frohs[plink.frohs$covg == c, 'bin4.froh'])/sqrt(15))
    OUT <- rbind(OUT, save)
  }
  OUT <- as.data.frame(OUT)
  colnames(OUT) <- c('covg','gt.mean','gt.se','pl.mean','pl.se','plink.mean','plink.se')
  
  ## 95% CIs
  lines(c(1:5)-offset, OUT$gt.mean, col = gt.col, lwd = lwd)
  points(c(1:5)-offset, OUT$gt.mean, col = gt.col, pch = 19)
  arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$gt.mean - OUT$gt.se*1.96),
         y1 = (OUT$gt.mean + OUT$gt.se*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5), OUT$pl.mean, col = pl.col, lwd = lwd)
  points(c(1:5), OUT$pl.mean, col = pl.col, pch = 19)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = (OUT$pl.mean - OUT$pl.se*1.96),
         y1 = (OUT$pl.mean + OUT$pl.se*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length=err.width)
  
  lines(c(1:5)+offset, OUT$plink.mean, col = plink.col, lwd = lwd)
  points(c(1:5)+offset, OUT$plink.mean, col = plink.col, pch = 19)
  arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$plink.mean - OUT$plink.se*1.96),
         y1 = (OUT$plink.mean + OUT$plink.se*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length=err.width)
  
  # legend('topright', pch = 19, lwd = lwd, legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'), bty = 'n', inset = c(0, 0.07), 
  #        col = c(gt.col, pl.col, plink.col), cex = txt.size)
  
  dev.off()
  
  k <- k+1
}

##### 3. Plotting f(ROH) by coverage, individual lines #####
alph <- 1  
bg.alph <- 0.25 ## for plot with mean lines +/- CIs
lwd <- 1.5
n <- 15
pt.size <- 1.5
err.wid <- 0.05
txt.size <- 1.25
ymin <- 0.075
ymax <- 0.425
xmin <- 0.85
xmax <- 5.15

k <- 1
while(k == 1){
pdf('../manuscript/r_scripts_AMH/figures_output/empirical/empirical_coverage_vs_fROH_indivlines.pdf', width = 10, height = 4)
par(mfrow=c(1,3))
plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
     xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - empirical', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n') 
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
  axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
  for(i in unique(froh.results$num.id)){
    temp <- froh.results[froh.results$num.id == i,]
    temp <- temp[order(temp$covg),]
  
    lines(c(1:5), c(temp$gt.froh), col = alpha(gt.col, alph))
    points(c(1:5), c(temp$gt.froh), pch = 19, col = gt.col)
  }  

plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
       xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - empirical', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n') 
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
  axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
  for(i in unique(froh.results$num.id)){
    temp <- froh.results[froh.results$num.id == i,]
    temp <- temp[order(temp$covg),]
    
    lines(c(1:5), c(temp$pl.froh), col = alpha(pl.col, alph))
    points(c(1:5), c(temp$pl.froh), pch = 19, col = pl.col)
  }  
  
plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
       xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - empirical', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
  axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
    for(i in unique(froh.results$num.id)){
    temp <- froh.results[froh.results$num.id == i,]
    temp <- temp[order(temp$covg),]
    
    lines(c(1:5), c(temp$plink.froh), col = alpha(plink.col, alph))
    points(c(1:5), c(temp$plink.froh), pch = 19, col = plink.col)
  }  

### background individual lines + means and 95% CIs
pdf('../manuscript/r_scripts_AMH/figures_output/empirical/empirical_coverage_vs_fROH_indivlines_95CIs.pdf', width = 10, height = 4)
par(mfrow=c(1,3))
OUT <- NULL
for(c in sort(unique(froh.results$covg))){
  save <- c(c, 
            mean(froh.results[froh.results$covg == c, 'pl.froh']), sd(froh.results[froh.results$covg == c, 'pl.froh'])/sqrt(n),
            mean(froh.results[froh.results$covg == c, 'gt.froh']), sd(froh.results[froh.results$covg == c, 'gt.froh'])/sqrt(n),
            mean(froh.results[froh.results$covg == c, 'plink.froh']), sd(froh.results[froh.results$covg == c, 'plink.froh'])/sqrt(n))
  OUT <- rbind(OUT, save)
}
OUT <- as.data.frame(OUT)
colnames(OUT) <- c('covg','mean.pl','se.pl','mean.gt','se.gt','mean.plink','se.plink')
write.csv(OUT, '../manuscript/r_scripts_AMH/tables_output/covg_fROH_means_SEs_empirical.csv')

## GT 
plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
     xaxt = 'n', xlab = 'Coverage', main = 'Genotypes', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
  axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
  for(i in unique(froh.results$num.id)){
    temp <- froh.results[froh.results$num.id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), c(temp$gt.froh), col = alpha(gt.col, bg.alph))
  }  
  
  lines(c(1:5), c(OUT$mean.gt), col = gt.col, lwd = lwd)
  points(c(1:5), OUT$mean.gt, pch = 19, col = gt.col, cex = pt.size)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = c(OUT$mean.gt - OUT$se.gt*1.96),
         y1 = c(OUT$mean.gt + OUT$se.gt*1.96),
         lwd = lwd, col = gt.col, code=3, angle=90, length= err.wid)

## PL 
plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
     xaxt = 'n', xlab = 'Coverage', main = 'Likelihoods', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
  axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
  for(i in unique(froh.results$num.id)){
    temp <- froh.results[froh.results$num.id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), c(temp$pl.froh), col = alpha(pl.col, bg.alph))
  }  
  lines(c(1:5), c(OUT$mean.pl), col = pl.col, lwd = lwd)
  points(c(1:5), OUT$mean.pl, pch = 19, col = pl.col, cex = pt.size)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = c(OUT$mean.pl - OUT$se.pl*1.96),
         y1 = c(OUT$mean.pl + OUT$se.pl*1.96),
         lwd = lwd, col = pl.col, code=3, angle=90, length= err.wid)

## PLINK 
plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
     xaxt = 'n', xlab = 'Coverage', main = 'PLINK', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
  axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
  axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
  for(i in unique(froh.results$num.id)){
    temp <- froh.results[froh.results$num.id == i,]
    temp <- temp[order(temp$covg),]
    lines(c(1:5), c(temp$plink.froh), col = alpha(plink.col, bg.alph))
  }  
  
  lines(c(1:5), c(OUT$mean.plink), col = plink.col, lwd = lwd)
  points(c(1:5), OUT$mean.plink, pch = 19, col = plink.col, cex = pt.size)
  arrows(x0 = c(1:5), x1 = c(1:5), y0 = c(OUT$mean.plink - OUT$se.plink*1.96),
         y1 = c(OUT$mean.plink + OUT$se.plink*1.96),
         lwd = lwd, col = plink.col, code=3, angle=90, length= err.wid)

dev.off()

### background individual lines + means and 83% CIs
pdf('../manuscript/r_scripts_AMH/figures_output/empirical/empirical_coverage_vs_fROH_indivlines_83CIs.pdf', width = 10, height = 4)
par(mfrow=c(1,3))
OUT <- NULL
for(c in sort(unique(froh.results$covg))){
  save <- c(c, 
            mean(froh.results[froh.results$covg == c, 'pl.froh']), sd(froh.results[froh.results$covg == c, 'pl.froh'])/sqrt(n),
            mean(froh.results[froh.results$covg == c, 'gt.froh']), sd(froh.results[froh.results$covg == c, 'gt.froh'])/sqrt(n),
            mean(froh.results[froh.results$covg == c, 'plink.froh']), sd(froh.results[froh.results$covg == c, 'plink.froh'])/sqrt(n))
  OUT <- rbind(OUT, save)
}
OUT <- as.data.frame(OUT)
colnames(OUT) <- c('covg','mean.pl','se.pl','mean.gt','se.gt','mean.plink','se.plink')

## GT 
plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
     xaxt = 'n', xlab = 'Coverage', main = 'Genotypes', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
for(i in unique(froh.results$num.id)){
  temp <- froh.results[froh.results$num.id == i,]
  temp <- temp[order(temp$covg),]
  lines(c(1:5), c(temp$gt.froh), col = alpha(gt.col, bg.alph))
}  

lines(c(1:5), c(OUT$mean.gt), col = gt.col, lwd = lwd)
points(c(1:5), OUT$mean.gt, pch = 19, col = gt.col, cex = pt.size)
arrows(x0 = c(1:5), x1 = c(1:5), y0 = c(OUT$mean.gt - OUT$se.gt*1.37),
       y1 = c(OUT$mean.gt + OUT$se.gt*1.37),
       lwd = lwd, col = gt.col, code=3, angle=90, length= err.wid)

## PL 
plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
     xaxt = 'n', xlab = 'Coverage', main = 'Likelihoods', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
for(i in unique(froh.results$num.id)){
  temp <- froh.results[froh.results$num.id == i,]
  temp <- temp[order(temp$covg),]
  lines(c(1:5), c(temp$pl.froh), col = alpha(pl.col, bg.alph))
}  
lines(c(1:5), c(OUT$mean.pl), col = pl.col, lwd = lwd)
points(c(1:5), OUT$mean.pl, pch = 19, col = pl.col, cex = pt.size)
arrows(x0 = c(1:5), x1 = c(1:5), y0 = c(OUT$mean.pl - OUT$se.pl*1.37),
       y1 = c(OUT$mean.pl + OUT$se.pl*1.37),
       lwd = lwd, col = pl.col, code=3, angle=90, length= err.wid)

## PLINK 
plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)', 
     xaxt = 'n', xlab = 'Coverage', main = 'PLINK', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0.0','0.1','0.2','0.3','0.4'), cex.axis = txt.size)
for(i in unique(froh.results$num.id)){
  temp <- froh.results[froh.results$num.id == i,]
  temp <- temp[order(temp$covg),]
  lines(c(1:5), c(temp$plink.froh), col = alpha(plink.col, bg.alph))
}  

lines(c(1:5), c(OUT$mean.plink), col = plink.col, lwd = lwd)
points(c(1:5), OUT$mean.plink, pch = 19, col = plink.col, cex = pt.size)
arrows(x0 = c(1:5), x1 = c(1:5), y0 = c(OUT$mean.plink - OUT$se.plink*1.37),
       y1 = c(OUT$mean.plink + OUT$se.plink*1.37),
       lwd = lwd, col = plink.col, code=3, angle=90, length= err.wid)

dev.off()

k <- k+1
}



##### 99. Plotting individual called ROHs for all 3 analyses and coverage levels #####
## set up data for plotting
gt.chrom.plots <- gt.res
pl.chrom.plots <- pl.res
plink.chrom.plots <- plink.res
## plot gap
gap <- 100e6
plot.tot.len <- tot.len + 5*gap

gt.chrom.plots$plot.start <- gt.chrom.plots$start
gt.chrom.plots$plot.end <- gt.chrom.plots$end
for(c in c(1:5)){
  gt.chrom.plots[gt.chrom.plots$chrom.num == c, 'plot.start'] <- gt.chrom.plots[gt.chrom.plots$chrom.num == c, 'start'] + chroms$cum.len[c-1] + gap*(c-1)
  gt.chrom.plots[gt.chrom.plots$chrom.num == c, 'plot.end'] <- gt.chrom.plots[gt.chrom.plots$chrom.num == c, 'end'] + chroms$cum.len[c-1] + gap*(c-1 )
}
pl.chrom.plots$plot.start <- pl.chrom.plots$start
pl.chrom.plots$plot.end <- pl.chrom.plots$end
for(c in c(1:5)){
  pl.chrom.plots[pl.chrom.plots$chrom.num == c, 'plot.start'] <- pl.chrom.plots[pl.chrom.plots$chrom.num == c, 'start'] + chroms$cum.len[c-1] + gap*(c-1)
  pl.chrom.plots[pl.chrom.plots$chrom.num == c, 'plot.end'] <- pl.chrom.plots[pl.chrom.plots$chrom.num == c, 'end'] + chroms$cum.len[c-1] + gap*(c-1 )
}
plink.chrom.plots$plot.start <- plink.chrom.plots$start
plink.chrom.plots$plot.end <- plink.chrom.plots$end
for(c in c(1:5)){
  plink.chrom.plots[plink.chrom.plots$chrom == c, 'plot.start'] <- plink.chrom.plots[plink.chrom.plots$chrom == c, 'start'] + chroms$cum.len[c-1] + gap*(c-1)
  plink.chrom.plots[plink.chrom.plots$chrom == c, 'plot.end'] <- plink.chrom.plots[plink.chrom.plots$chrom == c, 'end'] + chroms$cum.len[c-1] + gap*(c-1 )
}

wid <- 3 ## line widths in plot ~ or ~
hit <- 0.1 ## polygon height

## plot all chromosomes
pdf('../manuscript/r_scripts_AMH/figures_output/empirical/individual_true_and_called_ROHs_across_all_chromosomes.pdf', width = 15, height = 7)
par(mar = c(5.1, 5.1, 4.1, 2.1))

for(i in unique(pl.chrom.plots$id)){
  sub.pl <- pl.chrom.plots[pl.chrom.plots$id == i,]
  sub.gt <- gt.chrom.plots[gt.chrom.plots$id == i,]
  sub.pk <- plink.chrom.plots[plink.chrom.plots$id == i,]
  
  plot(0,0, xlim = c(1, plot.tot.len), ylim = c(0, 15), col = 'transparent', yaxt = 'n', xlab = '', ylab = '', main = i, xaxt = 'n', bty = 'n')
    lines(x = c(-1e9, 4e9), y = c(5.5, 5.5), lty = 2)
    lines(x = c(-1e9, 4e9), y = c(10.5, 10.5), lty = 2)
    axis(2, at = c(1:15), labels = c('5X','10X','15X','30X','50X','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X'), las = 2)  
    mtext('PLINK', side = 2, line = 3, at = 13)
    mtext('Genotypes only', side = 2, line = 3, at = 8)
    mtext('Genotype\nlikelihoods', side = 2, line = 3, at = 3)
    for(c in c(1:6)){
      lines(x = c((chroms$cum.len[c] + (gap * (c-1)) - chroms$length[c] + 1), (chroms$cum.len[c] + (gap * (c-1)))), y = c(0.5,0.5))
    }
    for(c in c(1:6)){
      text(x = ((chroms$cum.len[c] + (gap * (c-1))) - chroms$length[c]/2), y = 0, labels = c)
    }
    
  y <- 1
  for(c in c(5, 10, 15, 30, 50)){
    temp <- sub.pl[sub.pl$covg == c,]
    if(nrow(temp) > 0){
      for(r in 1:nrow(temp)){
        polygon(x = c(temp$plot.start[r], temp$plot.end[r], temp$plot.end[r], temp$plot.start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = pl.col, border = NA)
      }
    }
    y <- y+1
  }
  for(c in c(5, 10, 15, 30, 50)){
    temp <- sub.gt[sub.gt$covg == c,]
    if(nrow(temp) > 0){
      for(r in 1:nrow(temp)){
        polygon(x = c(temp$plot.start[r], temp$plot.end[r], temp$plot.end[r], temp$plot.start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = gt.col, border = NA)
      }
    }
    y <- y+1
  }
  for(c in c(5, 10, 15, 30, 50)){
    temp <- sub.pk[sub.pk$covg == c,]
    if(nrow(temp) > 0){
      for(r in 1:nrow(temp)){
        polygon(x = c(temp$plot.start[r], temp$plot.end[r], temp$plot.end[r], temp$plot.start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = plink.col, border = NA)
      }
    }
    y <- y+1
  }
}
dev.off()

## just plot first 10 Mb of chromosome 1
pdf('../manuscript/r_scripts_AMH/figures_output/empirical/individual_true_and_called_ROHs_across_chromosome_1_10Mb_window.pdf', width = 15, height = 7)
par(mar = c(5.1, 5.1, 4.1, 2.1))

for(i in unique(pl.chrom.plots$id)){
  sub.pl <- pl.chrom.plots[pl.chrom.plots$id == i,]
  sub.gt <- gt.chrom.plots[gt.chrom.plots$id == i,]
  sub.pk <- plink.chrom.plots[plink.chrom.plots$id == i,]
  
  plot(0,0, xlim = c(1, 10e6), ylim = c(0, 15), col = 'transparent', yaxt = 'n', xlab = '', ylab = '', main = i, xaxt = 'n', bty = 'n')
  lines(x = c(-1e9, 4e9), y = c(5.5, 5.5), lty = 2)
  lines(x = c(-1e9, 4e9), y = c(10.5, 10.5), lty = 2)
  axis(2, at = c(1:15), labels = c('5X','10X','15X','30X','50X','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X'), las = 2)  
  mtext('PLINK', side = 2, line = 3, at = 13)
  mtext('Genotypes only', side = 2, line = 3, at = 8)
  mtext('Genotype\nlikelihoods', side = 2, line = 3, at = 3)
  for(c in c(1:6)){
    lines(x = c((chroms$cum.len[c] + (gap * (c-1)) - chroms$length[c] + 1), (chroms$cum.len[c] + (gap * (c-1)))), y = c(0.5,0.5))
  }
  for(c in c(1:6)){
    text(x = ((chroms$cum.len[c] + (gap * (c-1))) - chroms$length[c]/2), y = 0, labels = c)
  }
  
  y <- 1
  for(c in c(5, 10, 15, 30, 50)){
    temp <- sub.pl[sub.pl$covg == c,]
    if(nrow(temp) > 0){
      for(r in 1:nrow(temp)){
        polygon(x = c(temp$plot.start[r], temp$plot.end[r], temp$plot.end[r], temp$plot.start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = pl.col, border = NA)
      }
    }
    y <- y+1
  }
  for(c in c(5, 10, 15, 30, 50)){
    temp <- sub.gt[sub.gt$covg == c,]
    if(nrow(temp) > 0){
      for(r in 1:nrow(temp)){
        polygon(x = c(temp$plot.start[r], temp$plot.end[r], temp$plot.end[r], temp$plot.start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = gt.col, border = NA)
      }
    }
    y <- y+1
  }
  for(c in c(5, 10, 15, 30, 50)){
    temp <- sub.pk[sub.pk$covg == c,]
    if(nrow(temp) > 0){
      for(r in 1:nrow(temp)){
        polygon(x = c(temp$plot.start[r], temp$plot.end[r], temp$plot.end[r], temp$plot.start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = plink.col, border = NA)
      }
    }
    y <- y+1
  }
}
dev.off()
