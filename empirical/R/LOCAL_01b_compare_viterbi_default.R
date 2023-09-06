library(scales)
library(ghibli)
library(TeachingDemos)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/data/')

### set colors
gt.col <- ghibli_palette('PonyoMedium')[4]
lt.gt.col <- ghibli_palette('PonyoLight')[4]
pl.col <- ghibli_palette('PonyoMedium')[2]
lt.pl.col <- ghibli_palette('PonyoLight')[2]
min.covg.col <- ghibli_palette('PonyoMedium')[6]
max.covg.col <- ghibli_palette('PonyoMedium')[1]

pl.pal <- colorRampPalette(c(lt.pl.col, max.covg.col))
pl.cols <- pl.pal(5)
gt.pal <- colorRampPalette(c(gt.col, max.covg.col))
gt.cols <- gt.pal(5)

## can set up defaults/vtrained colors instead of focusing on coverage-based shades
pl.def.col <- pl.cols[3]
pl.vit.col <- pl.cols[1]
gt.def.col <- gt.cols[3]
gt.vit.col <- gt.cols[1]

`%notin%` <- Negate(`%in%`)

## tas dev chrom data
chroms <- read.table('mSarHar1.11_autosomes.txt', sep = '\t', header = TRUE)
chroms <- chroms[,c(5,9)]
colnames(chroms) <- c('chrom','length')
tot.len <- sum(chroms$length)

## set empirical ROH length bin cutoffs and names
b.1 <- 100e3
b.2 <- 500e3
b.3 <- 1e6
b.4 <- 2e6

bin.names <- c('Short','Intermediate','Long','Very long')

## summarize individual f(ROH) values (overall and length-specific) for comparisons
OUT <- NULL
for(m in c('GT','PL')){
  dat <- read.csv(paste0('tasdev_',m,'_allcovg_levels.csv'))
  for(i in unique(dat$id)){
    sub <- dat[dat$id == i,]
    for(c in unique(sub$covg)){
      for(e in unique(sub$ests)){
        froh <- sum(sub[sub$covg == c & sub$ests == e, 'length'])/tot.len
        froh.1 <- sum(sub[sub$covg == c & sub$ests == e &
                          sub$length >= b.1 & sub$length < b.2, 'length'])/tot.len
        froh.2 <- sum(sub[sub$covg == c & sub$ests == e &
                            sub$length >= b.2 & sub$length < b.3, 'length'])/tot.len
        froh.3 <- sum(sub[sub$covg == c & sub$ests == e &
                            sub$length >= b.3 & sub$length < b.4, 'length'])/tot.len
        froh.4 <- sum(sub[sub$covg == c & sub$ests == e &
                            sub$length >= b.4, 'length'])/tot.len
        save <- c(m, i, c, e, froh, froh.1, froh.2, froh.3, froh.4)
        OUT <- rbind(OUT, save)
      }
    }
  }
}

dat <- as.data.frame(OUT)
for(c in c(5:9)){
  dat[,c] <- as.numeric(dat[,c])
}
colnames(dat) <- c('method','id','covg','hmm','froh','froh.1','froh.2','froh.3','froh.4')

for(m in unique(dat$method)){
  sub <- dat[dat$method == m,]
  cats <- expand.grid(unique(sub$ests), unique(sub$covg))
  
  ## plot how changing HMM transition probabilites affects overall f(ROH)
  pdf(paste0('../figures/',m,'_viterbi_vs_default_fROH.pdf'), width = 7, height = 5)
  plot(0,0, main = m,
       xlim = c(1, c(length(unique(sub$covg)) * length(unique(sub$ests)))), xaxt = 'n', xlab = '',
       ylim = c(0, max(dat$froh)), ylab = 'f(ROH)',
       col = 'transparent')
    axis(1, at = c(1:c(length(unique(sub$covg)) * length(unique(sub$ests)))),
         labels = paste0(cats[,2],'X\n',cats[,1]), tick = FALSE)
    x <- 1
    for(c in unique(sub$covg)){
      for(i in unique(sub$id)){
        lines(c(x, x+1), 
              c(sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[1], 'froh'],
                sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[2], 'froh']),
              col = alpha('darkgrey', 0.6))
        points(c(x, x+1), 
              c(sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[1], 'froh'],
                sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[2], 'froh']),
              col = alpha('black', 1), pch = 16, cex = 0.75)
      }
      x <- x+2
    }
    dev.off()
    
    ## " length-specific f(ROH)
    pdf(paste0('../figures/',m,'_viterbi_vs_default_fROH_by_bins.pdf'), width = 7, height = 5)
    n <- 1
    for(l in c(6:9)){
      plot(0,0, main = bin.names[n],
           xlim = c(1, c(length(unique(sub$covg)) * length(unique(sub$ests)))), xaxt = 'n', xlab = '',
           ylim = c(0, max(dat[,c(6:9)])), ylab = 'f(ROH)',
           col = 'transparent')
      axis(1, at = c(1:c(length(unique(sub$covg)) * length(unique(sub$ests)))),
           labels = paste0(cats[,2],'X\n',cats[,1]), tick = FALSE)
      x <- 1
      for(c in unique(sub$covg)){
        for(i in unique(sub$id)){
          lines(c(x, x+1), 
                c(sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[1], l],
                  sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[2], l]),
                col = alpha('darkgrey', 0.6))
          points(c(x, x+1), 
                 c(sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[1], l],
                   sub[sub$id == i & sub$covg == c & sub$ests == unique(sub$ests)[2], l]),
                 col = alpha('black', 1), pch = 16, cex = 0.75)
        }
        x <- x+2
      }
      n <- n+1
    }
    dev.off()
}

##### Mirror Fig. 6 but for Viterbi vs. defaults (line plot per method x ROH length bin) #####
xmin <- 0.85
xmax <- 5.15
ymin <- 0
ymax <- 0.3
txt.size <- 1.25
bg.alph <- 0.2
offset <- 0.05
lwd <- 2
pt.cex <- 1
err.width <- 0.1

variables <- as.data.frame(cbind(c('GT','PL'),
                                 c(gt.def.col, pl.def.col),
                                 c(gt.vit.col, pl.vit.col)))
colnames(variables) <- c('method','def.col','vit.col')

bins <- cbind(bin.names, c(1:4), c(6:9))

for(r in 1:nrow(variables)){
  sub.dat <- dat[dat$method == variables$method[r],]
  sub.dat$covg <- as.numeric(sub.dat$covg)
  def.col <- variables$def.col[r]
  vit.col <- variables$vit.col[r]
  method <- variables$method[r]
  
  pdf(paste0('../figures/',method,'_viterbi_vs_default_fROH_by_length_bins_indivlines_95_CIs.pdf'), width = 4.75, height = 6)
  par(mfrow = c(1,1), mar = c(5.1, 4.6, 4.1, 2.1))
  
  for(b in 1:nrow(bins)){
    colm <- as.numeric(bins[b,3])
    
    plot(0,0, xlim = c(xmin,xmax), ylim = c(ymin, ymax), 
         xaxt = 'n', main = bins[b,1], xlab = 'Coverage', ylab = substitute(paste(italic('F')[ROH])),
         cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
    axis(1, at = c(1,2,3,4,5), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
    axis(2, at = seq(0, ymax, 0.05), labels = c('0.0','','0.1','','0.2','','0.3'), cex.axis = txt.size)
    x <- 1
    
    for(i in unique(dat$id)){
      ## default
      temp <- sub.dat[sub.dat$id == i & sub.dat$hmm == 'default',]
      temp <- temp[order(temp$covg),]
      lines(c(1:5), temp[,colm], col = alpha(def.col, bg.alph))
      ## viterbi
      temp <- sub.dat[sub.dat$id == i & sub.dat$hmm == 'vtrained',]
      temp <- temp[order(temp$covg),]
      lines(c(1:5), temp[,colm], col = alpha(vit.col, bg.alph))
    }
    OUT <- NULL
    for(c in sort(unique(sub.dat$covg))){
      save <- c(c, 
                mean(sub.dat[sub.dat$covg == c & sub.dat$hmm == 'default', colm]), 
                sd(sub.dat[sub.dat$covg == c & sub.dat$hmm == 'default', colm])/sqrt(15),
                mean(sub.dat[sub.dat$covg == c & sub.dat$hmm == 'vtrained', colm]), 
                sd(sub.dat[sub.dat$covg == c & sub.dat$hmm == 'vtrained', colm])/sqrt(15))
      OUT <- rbind(OUT, save)
    }
    OUT <- as.data.frame(OUT)
    colnames(OUT) <- c('covg','def.mean','def.se','vit.mean','vit.se')
    write.csv(OUT, paste0('../output/',method,'_empirical_bin',bins[b,2],'.csv'), row.names = FALSE)
    
    ## 95% CIs
    lines(c(1:5)-offset, OUT$def.mean, col = def.col, lwd = lwd)
    points(c(1:5)-offset, OUT$def.mean, col = def.col, pch = 16, cex = pt.cex)
    suppressWarnings(arrows(x0 = c(1:5)-offset, x1 = c(1:5)-offset, y0 = (OUT$def.mean - OUT$def.se*1.96),
           y1 = (OUT$def.mean + OUT$def.se*1.96),
           lwd = lwd, col = def.col, code=3, angle=90, length=err.width))
    
    lines(c(1:5)+offset, OUT$vit.mean, col = vit.col, lwd = lwd)
    points(c(1:5)+offset, OUT$vit.mean, col = vit.col, pch = 16, cex = pt.cex)
    suppressWarnings(arrows(x0 = c(1:5)+offset, x1 = c(1:5)+offset, y0 = (OUT$vit.mean - OUT$vit.se*1.96),
           y1 = (OUT$vit.mean + OUT$vit.se*1.96),
           lwd = lwd, col = vit.col, code=3, angle=90, length=err.width))
    
    if(b == 4){
      legend('topright', legend = c('Default','Viterbi-trained'), col = c(def.col, vit.col), pch = 16, lwd = lwd,
             inset = 0.02, bty = 'n', cex = txt.size)
    }
  }
  dev.off()
}

