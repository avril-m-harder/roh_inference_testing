library(scales)
library(ghibli)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/')

### set colors
gt.col <- ghibli_palette('PonyoMedium')[4]
lt.gt.col <- ghibli_palette('PonyoLight')[4]
pl.col <- ghibli_palette('PonyoMedium')[2]
lt.pl.col <- ghibli_palette('PonyoLight')[2]
min.covg.col <- ghibli_palette('PonyoMedium')[6]
max.covg.col <- ghibli_palette('PonyoMedium')[1]

`%notin%` <- Negate(`%in%`)

##### Read in all ROH-calling results #####
res <- read.table('bcftools_output/bcftoolsROH_all_coordinates.txt', header = TRUE)

##### Loop over demographic scenarios #####
true.fns <- list.files(path = 'slim_true_data/', pattern = 'roh_coords')

## set ROH length bin cutoffs
b.1 <- 100e3
b.2 <- 500e3
b.3 <- 1e6
b.4 <- 2e6

OUT <- NULL

### remove when final data sets are in place:
true.fns <- true.fns[c(1,2,4,5)]
for(t in true.fns){
  demo <- unlist(strsplit(t, split = '_'))[1]
  print(demo)
  demo.res <- res[res$demo == demo,]
  demo.res$called.roh.id <- c(1:nrow(demo.res))
  overlap <- read.table(paste0('bcftools_output/bcftoolsROH_',demo,'_overlap_results.txt'), header = TRUE) 
  
  ## Read in known/true heterozygosity + ROH information
  true.rohs <- read.table(paste0('slim_true_data/',t))
  colnames(true.rohs) <- c('id','start','end','length')
  true.rohs <- true.rohs[true.rohs$length >= 100000,] ### only keeping true ROHs >= 100 kb because that's all we're evaluating the ability to call
  ## create unique ROH ID for linking true ROHs to called ROHs
  true.rohs$true.roh.id <- c(1:nrow(true.rohs))
  chrom.len <- 30e6
  
  for(m in unique(overlap$method)){
    for(c in unique(overlap$covg)){
      dat <- overlap[overlap$method == m & overlap$covg == c,]
      
      ## for each individual, calculate true and called f(ROH) values, ratio of true:called ROHs (i.e.,
      ## the ROH-lumping issue), and true and called length-specific f(ROH) values
      for(i in unique(dat$id)){
        sub.def <- dat[dat$id == i & dat$hmm == 'defaults',]
        sub.vit <- dat[dat$id == i & dat$hmm == 'vtrained',]
        sub.true <- true.rohs[true.rohs$id == i,]
        
        def.true.call.ratio <- length(unique(sub.def$true.roh.id))/length(unique(sub.def$called.roh.id))
        vit.true.call.ratio <- length(unique(sub.vit$true.roh.id))/length(unique(sub.vit$called.roh.id))
        
        uniq.def <- sub.def[!duplicated(sub.def$called.roh.id),]
        uniq.vit <- sub.vit[!duplicated(sub.vit$called.roh.id),]
        
        d.froh <- sum(uniq.def[!duplicated(uniq.def$called.roh.id), 'len'])/chrom.len
        v.froh <- sum(uniq.vit[!duplicated(uniq.vit$called.roh.id), 'len'])/chrom.len
        true.froh <- sum(sub.true$length)/chrom.len
        
        d.bin1 <- sum(uniq.def[uniq.def$len >= b.1 & uniq.def$len < b.2, 'len'])/chrom.len
        d.bin2 <- sum(uniq.def[uniq.def$len >= b.2 & uniq.def$len < b.3, 'len'])/chrom.len
        d.bin3 <- sum(uniq.def[uniq.def$len >= b.3 & uniq.def$len < b.4, 'len'])/chrom.len
        d.bin4 <- sum(uniq.def[uniq.def$len >= b.4, 'len'])/chrom.len

        v.bin1 <- sum(uniq.vit[uniq.vit$len >= b.1 & uniq.vit$len < b.2, 'len'])/chrom.len
        v.bin2 <- sum(uniq.vit[uniq.vit$len >= b.2 & uniq.vit$len < b.3, 'len'])/chrom.len
        v.bin3 <- sum(uniq.vit[uniq.vit$len >= b.3 & uniq.vit$len < b.4, 'len'])/chrom.len
        v.bin4 <- sum(uniq.vit[uniq.vit$len >= b.4, 'len'])/chrom.len

        true.bin1 <- sum(sub.true[sub.true$length >= b.1 & sub.true < b.2, 'length'])/chrom.len
        true.bin2 <- sum(sub.true[sub.true$length >= b.2 & sub.true < b.3, 'length'])/chrom.len
        true.bin3 <- sum(sub.true[sub.true$length >= b.3 & sub.true < b.4, 'length'])/chrom.len
        true.bin4 <- sum(sub.true[sub.true$length >= b.4, 'length'])/chrom.len

        save <- c(demo, m, c, i, d.froh, v.froh, true.froh,
                  d.bin1, d.bin2, d.bin3, d.bin4,
                  v.bin1, v.bin2, v.bin3, v.bin4,
                  true.bin1, true.bin2, true.bin3, true.bin4,
                  def.true.call.ratio, vit.true.call.ratio)
        save[is.na(save)] <- 0
        OUT <- rbind(OUT, save)
      }
    }
  }
}

dat <- as.data.frame(OUT)
colnames(dat) <- c('demo','method', 'covg', 'id', 'd.froh', 'v.froh', 'true.froh',
                   'd.bin1', 'd.bin2', 'd.bin3', 'd.bin4',
                   'v.bin1', 'v.bin2', 'v.bin3', 'v.bin4',
                   'true.bin1', 'true.bin2', 'true.bin3', 'true.bin4',
                   'def.true.call.ratio','vit.true.call.ratio')
for(c in 5:ncol(dat)){
  dat[,c] <- as.numeric(dat[,c])
}

## plot results
pt.alph <- 0.9

OUT <- NULL ## place to store correlation coefficients

pdf('/Users/Avril/Desktop/test.pdf', width = 16, height = 4.5)
par(mfrow = c(1, 3))
for(d in unique(dat$demo)){
  for(m in unique(dat$method)){
    for(c in unique(dat$covg)){
      ## set color
      if(m == 'GT'){
        colour <- gt.col
        lt.colour <- lt.gt.col
      } else{
        colour <- pl.col
        lt.colour <- lt.pl.col
      }
      sub <- dat[dat$demo == d & dat$method == m & dat$covg == c,]
      
      def <- cor(sub$true.froh, sub$d.froh)
      vit <- cor(sub$true.froh, sub$v.froh)
      if(def > vit){
        better <- 'default'
      } else{
        better <- 'vtrained'
      }
      save <- c(d, m, c, def, vit, better, 
                mean(abs(sub$true.froh - sub$d.froh), na.rm = TRUE), sd(abs(sub$true.froh - sub$d.froh), na.rm = TRUE),
                mean(abs(sub$true.froh - sub$v.froh), na.rm = TRUE), sd(abs(sub$true.froh - sub$v.froh), na.rm = TRUE),
                mean(sub$def.true.call.ratio, na.rm = TRUE), sd(sub$def.true.call.ratio, na.rm = TRUE), 
                mean(sub$vit.true.call.ratio, na.rm = TRUE), sd(sub$vit.true.call.ratio, na.rm = TRUE))
      OUT <- rbind(OUT, save)
      
      ## true f(ROH) vs. default and viterbi f(ROH)s
      par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = FALSE)
      plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), 
           pch = 16, col = 'transparent',
           main = paste0(d,' - ',m,' - ',c), xlab = 'True f(ROH)', ylab = 'Called f(ROH)')
        abline(0, 1, lty = 2, col = 'darkgrey')
        points(sub$true.froh, sub$d.froh, pch = 16, col = alpha(colour, pt.alph))
        abline(lm(sub$d.froh ~ sub$true.froh), col = colour)
        points(sub$true.froh, sub$v.froh, pch = 17, col = alpha(lt.colour, pt.alph))
        abline(lm(sub$v.froh ~ sub$true.froh), col = lt.colour)
        par(xpd = TRUE)
        legend('right', pch = c(16, 17), col = c(colour, lt.colour), legend = c('Default','Viterbi-trained'), inset = c(-0.27, 0))
        
      par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = FALSE)
      plot(0, 0, xlim = c(min(sub$true.froh), max(sub$true.froh)),
           ylim = c(min(sub$d.froh, sub$v.froh), max(sub$d.froh, sub$v.froh)), 
           pch = 16, col = 'transparent',
           main = paste0(d,' - ',m,' - ',c), xlab = 'True f(ROH)', ylab = 'Called f(ROH)')
        abline(0, 1, lty = 2, col = 'darkgrey')
        points(sub$true.froh, sub$d.froh, pch = 16, col = alpha(colour, pt.alph))
        abline(lm(sub$d.froh ~ sub$true.froh), col = colour)
        points(sub$true.froh, sub$v.froh, pch = 17, col = alpha(lt.colour, pt.alph))
        abline(lm(sub$v.froh ~ sub$true.froh), col = lt.colour)
        par(xpd = TRUE)
        legend('right', pch = c(16, 17), col = c(colour, lt.colour), legend = c('Default','Viterbi-trained'), inset = c(-0.27, 0))
    
      # ## true binned f(ROH) vs. default binned f(ROH)
      # par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = FALSE)
      # plot(0, 0, pch = 16, col = 'transparent', xlim = c(0, 1), ylim = c(0, 1),
      #      main = paste0(d,' - ',m,' - ',c,' - Default'), xlab = 'True binned f(ROH)', ylab = 'Called binned f(ROH)')
      #   abline(0, 1, lty = 2, col = 'darkgrey')
      #   points(sub$true.bin1, sub$v.bin1, pch = 16, col = alpha(colour, pt.alph))
      #   points(sub$true.bin2, sub$d.bin2, pch = 17, col = alpha(colour, pt.alph))
      #   points(sub$true.bin3, sub$d.bin3, pch = 15, col = alpha(colour, pt.alph))
      #   points(sub$true.bin4, sub$d.bin4, pch = 18, col = alpha(colour, pt.alph))
      #   par(xpd = TRUE)
      #   legend('right', pch = c(16, 17, 15, 18), legend = c('Short','Intermediate','Long','Very long'), inset = c(-0.27, 1),
      #          col = colour)
      #   
      # ## true binned f(ROH) vs. viterbi binned f(ROH)
      # par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = FALSE)
      # plot(0, 0, pch = 16, col = 'transparent', xlim = c(0, 1), ylim = c(0, 1),
      #      main = paste0(d,' - ',m,' - ',c,' - Viterbi-trained'), xlab = 'True binned f(ROH)', ylab = 'Called binned f(ROH)')
      #   abline(0, 1, lty = 2, col = 'darkgrey')
      #   points(sub$true.bin1, sub$v.bin1, pch = 16, col = alpha(colour, pt.alph))
      #   points(sub$true.bin2, sub$v.bin2, pch = 17, col = alpha(colour, pt.alph))
      #   points(sub$true.bin3, sub$v.bin3, pch = 15, col = alpha(colour, pt.alph))
      #   points(sub$true.bin4, sub$v.bin4, pch = 18, col = alpha(colour, pt.alph))
      #   par(xpd = TRUE)
      #   legend('right', pch = c(16, 17, 15, 18), legend = c('Short','Intermediate','Long','Very long'), inset = c(-0.27, 1),
      #          col = colour)
        
      ## # true : # called ROHs
      
    }
  }
}
dev.off()

OUT <- OUT[,-6]
colnames(OUT) <- c('demo','method','covg','default.cor','viterbi.cor','mean.d.froh.diff','sd.d.froh.diff',
                   'mean.v.froh.diff','sd.v.froh.diff','mean.def.true.call.ratio','sd.def.true.call.ratio',
                   'mean.vit.true.call.ratio','sd.vit.true.call.ratio')

dat <- as.data.frame(OUT)
for(c in 4:ncol(dat)){
  dat[,c] <- as.numeric(dat[,c])
}

pdf('../figures/bcftools_default_vs_viterbi_froh_diff.pdf', width = 6, height = 6)
for(d in unique(dat$demo)){
  sub <- dat[dat$demo == d,]
  par(xpd = FALSE)
  plot(0, 0, col = 'transparent', xlim = c(1, 20), 
       # ylim = c(0, max(c(sub$mean.d.froh.diff+sub$sd.d.froh.diff, sub$mean.v.froh.diff+sub$sd.v.froh.diff))), 
       ylim = c(0, 0.7),
       xaxt = 'n', ylab = '', xlab = '', main = d)
    abline(v = 10.5, lty = 2, col = 'darkgrey')
    points(seq(1, 9, 2), sub[sub$method == 'GT', 'mean.d.froh.diff'], pch = 16, col = gt.col)
    for(c in 1:5){
      lines(c((c*2-1), (c*2-1)), 
            c((sub[sub$method == 'GT', 'mean.d.froh.diff'][c] - sub[sub$method == 'GT', 'sd.d.froh.diff'][c]),
            (sub[sub$method == 'GT', 'mean.d.froh.diff'][c] + sub[sub$method == 'GT', 'sd.d.froh.diff'][c])),
            col = gt.col, lwd = 2)
    }
    points(seq(2, 10, 2), sub[sub$method == 'GT', 'mean.v.froh.diff'], pch = 17, col = gt.col)
    for(c in 1:5){
      lines(c((c*2), (c*2)), 
            c((sub[sub$method == 'GT', 'mean.v.froh.diff'][c] - sub[sub$method == 'GT', 'sd.v.froh.diff'][c]),
              (sub[sub$method == 'GT', 'mean.v.froh.diff'][c] + sub[sub$method == 'GT', 'sd.v.froh.diff'][c])),
            col = gt.col, lwd = 2)
    }
    points(seq(11, 19, 2), sub[sub$method == 'PL', 'mean.d.froh.diff'], pch = 16, col = pl.col)
    for(c in 6:10){
      lines(c((c*2-1), (c*2-1)), 
            c((sub[sub$method == 'PL', 'mean.d.froh.diff'][c-5] - sub[sub$method == 'PL', 'sd.d.froh.diff'][c-5]),
              (sub[sub$method == 'PL', 'mean.d.froh.diff'][c-5] + sub[sub$method == 'PL', 'sd.d.froh.diff'][c-5])),
            col = pl.col, lwd = 2)
    }
    points(seq(12, 20, 2), sub[sub$method == 'PL', 'mean.v.froh.diff'], pch = 17, col = pl.col)
    for(c in 6:10){
      lines(c((c*2), (c*2)), 
            c((sub[sub$method == 'PL', 'mean.v.froh.diff'][c-5] - sub[sub$method == 'PL', 'sd.v.froh.diff'][c-5]),
              (sub[sub$method == 'PL', 'mean.v.froh.diff'][c-5] + sub[sub$method == 'PL', 'sd.v.froh.diff'][c-5])),
            col = pl.col, lwd = 2)
    }
    axis(1, at = seq(1.5, 19.5, 2), labels = c('5X','10X','15X','30X','50X','5X','10X','15X','30X','50X'))
    par(xpd = TRUE)
    legend('topleft', pch = c(16, 17), legend = c('Default','Viterbi-trained'), inset = c(0, -0.15), bty = 'n')
}
dev.off()

pdf('../figures/bcftools_default_vs_viterbi_truevcalled_ratio.pdf', width = 6, height = 6)
for(d in unique(dat$demo)){
  sub <- dat[dat$demo == d,]
  par(xpd = FALSE)
  plot(0, 0, col = 'transparent', xlim = c(1, 20), 
       # ylim = c(0, max(c(sub$mean.def.true.call.ratio+sub$sd.def.true.call.ratio, sub$mean.vit.true.call.ratio+sub$sd.vit.true.call.ratio))),
       ylim = c(0, 5.2),
       xaxt = 'n', ylab = '', xlab = '', main = d)
  abline(v = 10.5, lty = 2, col = 'darkgrey')
  points(seq(1, 9, 2), sub[sub$method == 'GT', 'mean.def.true.call.ratio'], pch = 16, col = gt.col)
  for(c in 1:5){
    lines(c((c*2-1), (c*2-1)), 
          c((sub[sub$method == 'GT', 'mean.def.true.call.ratio'][c] - sub[sub$method == 'GT', 'sd.def.true.call.ratio'][c]),
            (sub[sub$method == 'GT', 'mean.def.true.call.ratio'][c] + sub[sub$method == 'GT', 'sd.def.true.call.ratio'][c])),
          col = gt.col, lwd = 2)
  }
  points(seq(2, 10, 2), sub[sub$method == 'GT', 'mean.vit.true.call.ratio'], pch = 17, col = gt.col)
  for(c in 1:5){
    lines(c((c*2), (c*2)), 
          c((sub[sub$method == 'GT', 'mean.vit.true.call.ratio'][c] - sub[sub$method == 'GT', 'sd.vit.true.call.ratio'][c]),
            (sub[sub$method == 'GT', 'mean.vit.true.call.ratio'][c] + sub[sub$method == 'GT', 'sd.vit.true.call.ratio'][c])),
          col = gt.col, lwd = 2)
  }
  points(seq(11, 19, 2), sub[sub$method == 'PL', 'mean.def.true.call.ratio'], pch = 16, col = pl.col)
  for(c in 6:10){
    lines(c((c*2-1), (c*2-1)), 
          c((sub[sub$method == 'PL', 'mean.def.true.call.ratio'][c-5] - sub[sub$method == 'PL', 'sd.def.true.call.ratio'][c-5]),
            (sub[sub$method == 'PL', 'mean.def.true.call.ratio'][c-5] + sub[sub$method == 'PL', 'sd.def.true.call.ratio'][c-5])),
          col = pl.col, lwd = 2)
  }
  points(seq(12, 20, 2), sub[sub$method == 'PL', 'mean.vit.true.call.ratio'], pch = 17, col = pl.col)
  for(c in 6:10){
    lines(c((c*2), (c*2)), 
          c((sub[sub$method == 'PL', 'mean.vit.true.call.ratio'][c-5] - sub[sub$method == 'PL', 'sd.vit.true.call.ratio'][c-5]),
            (sub[sub$method == 'PL', 'mean.vit.true.call.ratio'][c-5] + sub[sub$method == 'PL', 'sd.vit.true.call.ratio'][c-5])),
          col = pl.col, lwd = 2)
  }
  axis(1, at = seq(1.5, 19.5, 2), labels = c('5X','10X','15X','30X','50X','5X','10X','15X','30X','50X'))
  par(xpd = TRUE)
  legend('topleft', pch = c(16, 17), legend = c('Default','Viterbi-trained'), inset = c(0, -0.15), bty = 'n')
}
dev.off()













