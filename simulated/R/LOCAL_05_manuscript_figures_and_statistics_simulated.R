##### Code for analyzing simulated data
library(scales)
library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggridges)
theme_set(theme_bw(16))
# source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")
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

### set scenario names
demo.names <- cbind(c('bottle','decline','large-2000','small','large-1000'),
                    c('Bottlenecked population','Declining population','Large population',
                      'Small population','Large population'))

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/')

##### Loop over demographic scenarios #####
true.fns <- list.files(path = 'slim_true_data/', pattern = 'true_roh_coords')
demos <- do.call(rbind, strsplit(true.fns, split = '_'))[,1]

for(d in demos){
  dir.create(paste0('../figures/',d), showWarnings = FALSE)
  scen.name <- demo.names[demo.names[,1] == d, 2]
  for(h in c('vtrained')){
    print(paste0(d,' - ',h))
    figure.ct <- 1
    ##### 1A. Read in data and summary files #####
    ### True ROH information
    true.rohs <- read.table(paste0('slim_true_data/',d,'_true_roh_coords.txt'))
    colnames(true.rohs) <- c('id','start','end','length','true.roh.id')
    true.rohs <- true.rohs[true.rohs$length >= 100000,]
    chrom.len <- 30e6
    if(d == 'decline'){
      true.rohs <- true.rohs[true.rohs$id %notin% c(33,48),]
    }

    ### Read in ROH calling results
    bcf.res <- read.table('bcftools_output/bcftoolsROH_all_coordinates.txt',
                          header = TRUE, sep = '\t')
    bcf.res$covg <- as.numeric(gsub('x', '', bcf.res$covg))
    bcf.res$id <- gsub('i', '', bcf.res$id)
    bcf.res <- bcf.res[bcf.res$demo == d,]
    if(d == 'decline'){
      bcf.res <- bcf.res[bcf.res$id %notin% c(33,48),]
    }

    ## GT
    bcf.gt.res <- bcf.res[bcf.res$method == 'GT' & bcf.res$hmm == h,]
    bcf.gt.res <- bcf.gt.res[bcf.gt.res$length >= 100000,] ## applying 100kb filter

    ## PL
    bcf.pl.res <- bcf.res[bcf.res$method == 'PL' & bcf.res$hmm == h,]
    bcf.pl.res <- bcf.pl.res[bcf.pl.res$length >= 100000,] ## applying 100kb filter

    ## PLINK
    plink.res <- read.table(paste0('plink_final_iteration/',d,'_PLINK_all_coordinates.txt'), header = TRUE)
    ##### !!! Update this for all demo scenarios, covgs if necessary !!! #####
    ## final selection = default settings
    plink.res <- plink.res[plink.res$phwh == 1 & plink.res$phwm == 5 & plink.res$phws == 50 &
                             plink.res$phwt == 0.05 & plink.res$phzs == 100 & plink.res$phzg == 1000,]
    plink.res$length <- plink.res$end - plink.res$start + 1
    plink.res$id <- gsub('i', '', plink.res$id)
    if(d == 'decline'){
      plink.res <- plink.res[plink.res$id %notin% c(33,48),]
    }
    plink.res$covg <- as.numeric(gsub('x', '', plink.res$covg))

    ### Write/read in called vs. true f(ROH) results
    OUT <- NULL
    for(i in unique(true.rohs$id)){
      true.froh <- sum(true.rohs[true.rohs$id == i, 'length'])/chrom.len
      for(c in unique(bcf.gt.res$covg)){
        gt.froh <- sum(bcf.gt.res[bcf.gt.res$id == i & bcf.gt.res$covg == c, 'length'])/chrom.len
        pl.froh <- sum(bcf.pl.res[bcf.pl.res$id == i & bcf.pl.res$covg == c, 'length'])/chrom.len
        pk.froh <- sum(plink.res[plink.res$id == i & plink.res$covg == c, 'length'])/chrom.len
        save <- c(i, c, true.froh, pl.froh, gt.froh, pk.froh)
        OUT <- rbind(OUT, save)
      }
    }
    colnames(OUT) <- c('id','covg','true.froh','pl.froh','gt.froh','plink.froh')
    write.csv(OUT, paste0('3_methods_results/',d,'_true_vs_called_froh_data.csv'), row.names = FALSE)
    froh.stats <- read.csv(paste0('3_methods_results/',d,'_true_vs_called_froh_data.csv'), header = TRUE)

    ### Called ROH data (i.e., overlap/true information - previously pl.out, gt.out, and plink.out in scripts, need to F+R)
    bcf.overlap <- read.table(paste0('bcftools_output/bcftoolsROH_',d,'_overlap_results.txt'), header = TRUE)
    if(d == 'decline'){
      bcf.overlap <- bcf.overlap[bcf.overlap$id %notin% c(33,48),]
    }
    bcf.overlap$covg <- as.numeric(gsub('x', '', bcf.overlap$covg))
    gt.overlap <- bcf.overlap[bcf.overlap$method == 'GT' & bcf.overlap$hmm == h,]
    pl.overlap <- bcf.overlap[bcf.overlap$method == 'PL' & bcf.overlap$hmm == h,]

    plink.overlap <- read.table(paste0('plink_final_iteration/',d,'_PLINK_overlap_results.txt'), sep = '\t', header = TRUE)
    plink.overlap$covg <- as.numeric(gsub('x', '', plink.overlap$covg))
    ##### !!! Update this for all demo scenarios, covgs if necessary !!! #####
    ## final selection = default settings
    plink.overlap <- plink.overlap[plink.overlap$phwh == 1 & plink.overlap$phwm == 5 & plink.overlap$phws == 50 & plink.overlap$phwt == 0.05 & plink.overlap$phzs == 100 & plink.overlap$phzg == 1000,]
    if(d == 'decline'){
      plink.overlap <- plink.overlap[plink.overlap$id %notin% c(33,48),]
    }

    ##### 1B. Plot histograms of ROH length distributions (true and called) #####
    alph <- 0.4
    x.max <- max(true.rohs$length, bcf.gt.res$length, bcf.pl.res$length, plink.res$length)

    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_ROH_length_histograms.pdf'), width = 6, height = 5)
    figure.ct <- figure.ct + 1

    ## true ROHs
    hist(true.rohs$length, breaks = 50, xlab = 'True ROH length', main = 'True simulated data', xlim = c(100e3, x.max))
    ## called BCFtools Genotypes
    hist(bcf.gt.res$length, breaks = 50, xlab = 'Called ROH length', main = 'BCFtools Genotypes - simulated data', xlim = c(100e3, x.max))
    ## called BCFtools Likelihoods
    hist(bcf.pl.res$length, breaks = 25, xlab = 'Called ROH length', main = 'BCFtools Likelihoods - simulated data', xlim = c(100e3, x.max))
    ## called PLINK
    hist(plink.res$length, breaks = 50, xlab = 'Called ROH length', main = 'PLINK - simulated data', xlim = c(100e3, x.max))

    plot(density(true.rohs$length), xlab = 'ROH length')
    # axis(1, at = c(0, 2.5e5, 5e5, 7.5e5, 1e6), labels = c(0, 250, 500, 750, 1000))
    polygon(density(plink.res$length), col = alpha(plink.col, alph), border = plink.col)
    polygon(density(bcf.pl.res$length), col = alpha(pl.col, alph), border = pl.col)
    polygon(density(bcf.gt.res$length), col = alpha(gt.col, alph), border = gt.col)
    legend('topright', legend = c('True','Genotypes', 'Likelihoods','PLINK'), border = c('black', gt.col, pl.col, plink.col),
           fill = c('white', alpha(gt.col, alph), alpha(pl.col, alph), alpha(plink.col, alph)), bty = 'n')

    plot(density(true.rohs$length), xlim = c(1e5, 4e6), xlab = 'ROH length')
    # axis(1, at = c(0, 1e5, 2e5, 3e5, 4e5, 5e5), labels = c(0, 100, 200, 300, 400, 500))
    polygon(density(plink.res$length), col = alpha(plink.col, alph), border = plink.col)
    polygon(density(bcf.pl.res$length), col = alpha(pl.col, alph), border = pl.col)
    polygon(density(bcf.gt.res$length), col = alpha(gt.col, alph), border = gt.col)
    legend('topright', legend = c('True','Genotypes', 'Likelihoods','PLINK'), border = c('black', gt.col, pl.col, plink.col),
           fill = c('white', alpha(gt.col, alph), alpha(pl.col, alph), alpha(plink.col, alph)), bty = 'n')
    dev.off()

    ##### 1C. Plot cumulative f(ROH) for each sample #####
    alph <- 0.4
    pt.cex <- 0.2
    maximum <- max(bcf.gt.res$length, bcf.pl.res$length, plink.res$length)
    y.max <- max(froh.stats$pl.froh, froh.stats$gt.froh, froh.stats$plink.froh)

    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_cumulative_fROH_figures_simulated.pdf'), width = 20, height = 12)
    figure.ct <- figure.ct + 1
    par(mfrow = c(3, 5))

    ## Genotypes
    for(c in unique(bcf.gt.res$covg)){
      plot(0, 0, col = 'transparent', xlim = c(1e5, maximum), ylim = c(0, y.max),
           main = paste0('Simulated: BCFtools Genotypes - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
      axis(1, at = seq(from = 0, to = maximum, by = 1e6), labels = seq(from = 0, to = maximum, by = 1e6)/1e6)
      for(i in unique(bcf.gt.res$id)){
        sub <- bcf.gt.res[bcf.gt.res$id == i & bcf.gt.res$covg == c,]
        sub <- sub[order(sub$length),]
        sub$cumlen <- cumsum(sub$length)
        sub$cumfroh <- sub$cumlen / chrom.len
        lines(sub$length, sub$cumfroh, col = alpha(gt.col, alph))
        points(sub$length, sub$cumfroh, col = gt.col, cex = pt.cex)
      }
    }

    ## Likelihoods
    for(c in unique(bcf.pl.res$covg)){
      plot(0, 0, col = 'transparent', xlim = c(1e5, maximum), ylim = c(0, y.max),
           main = paste0('Simulated: BCFtools Likelihoods - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
      axis(1, at = seq(from = 0, to = maximum, by = 1e6), labels = seq(from = 0, to = maximum, by = 1e6)/1e6)
      for(i in unique(bcf.pl.res$id)){
        sub <- bcf.pl.res[bcf.pl.res$id == i & bcf.pl.res$covg == c,]
        sub <- sub[order(sub$length),]
        sub$cumlen <- cumsum(sub$length)
        sub$cumfroh <- sub$cumlen / chrom.len
        lines(sub$length, sub$cumfroh, col = alpha(pl.col, alph))
        points(sub$length, sub$cumfroh, col = pl.col, cex = pt.cex)
      }
    }

    ## PLINK
    for(c in unique(plink.res$covg)){
      plot(0, 0, col = 'transparent', xlim = c(1e5, maximum), ylim = c(0, y.max),
           main = paste0('Simulated: PLINK - ', c, 'X'), xlab = 'ROH length (Mb)', ylab = 'Cumulative f(ROH)', xaxt = 'n')
      axis(1, at = seq(from = 0, to = maximum, by = 1e6), labels = seq(from = 0, to = maximum, by = 1e6)/1e6)
      for(i in unique(plink.res$id)){
        sub <- plink.res[plink.res$id == i & plink.res$covg == c,]
        sub <- sub[order(sub$length),]
        sub$cumlen <- cumsum(sub$length)
        sub$cumfroh <- sub$cumlen / chrom.len
        lines(sub$length, sub$cumfroh, col = alpha(plink.col, alph))
        points(sub$length, sub$cumfroh, col = plink.col, cex = pt.cex)
      }
    }

    dev.off()


    ##### 3. Calculation of false-positive and false-negative rates #####
    ##### >>> 3A. BCFtools/ROH - PL (genotype likelihoods) #####
    OUT <- NULL
    for(i in unique(bcf.pl.res$id)){
      print(i)
      for(c in unique(bcf.pl.res$covg)){
        ## get all called ROHs >= 100 kb
        called <- bcf.pl.res[bcf.pl.res$id == i & bcf.pl.res$covg == c,]
        tot.call.len <- sum(called[!duplicated(called$called.roh.id), 'length'])
        n.called <- length(unique(called$called.roh.id))

        ## get all true ROHs >= 100 kb
        true <- true.rohs[true.rohs$id == i & true.rohs$length >= 100000,]
        tot.true.len <- sum(true$length)
        n.true <- length(unique(true$true.roh.id))

        ## for each called ROH,
        TOT.TRUE <- 0
        for(r in unique(called$called.roh.id)){
          ## total up length that overlaps any >= 100 kb true ROH == tp == true positive
          true.call.len <- sum(pl.overlap[pl.overlap$called.roh.id == r, 'true.len'])
          TOT.TRUE <- TOT.TRUE + true.call.len
        }
        save <- c(i, c, tot.call.len, n.called, tot.true.len, n.true, TOT.TRUE)
        OUT <- rbind(OUT, as.integer(save))
      }
    }
    pl.true.v.called <- as.data.frame(OUT)
    colnames(pl.true.v.called) <- c('id','covg','total.call.len','n.called','tot.true.len','n.true','tp')

    ## total.call.len = total length of all called ROHs
    ## n.called = # of ROHs called
    ## tot.true.len = total length of all true ROHs
    ## n.true = # of true ROHs
    ## tp = # of called ROH base pairs that overlap true ROH base pairs

    for(c in 1:ncol(pl.true.v.called)){
      pl.true.v.called[,c] <- as.integer(pl.true.v.called[,c])
    }
    ## fn = false negative = total length of true ROHs - total length of called ROHs that overlap true
    ## ROHs = total length of true ROHs - true positive calls
    pl.true.v.called$fn <- pl.true.v.called$tot.true.len - pl.true.v.called$tp
    ## fp = false positive = total length of called ROHs - total length of called ROHs that overlap true
    ## ROHs = total length of called ROHs - true positive calls
    pl.true.v.called$fp <- pl.true.v.called$total.call.len - pl.true.v.called$tp
    ## tn = true negative = total chromosome length not covered by true positives, false negatives, or false positives
    pl.true.v.called$tn <- chrom.len - pl.true.v.called$tp - pl.true.v.called$fn - pl.true.v.called$fp
    pl.true.v.called$false.pos.rate <- pl.true.v.called$fp / (pl.true.v.called$fp + pl.true.v.called$tn)
    pl.true.v.called$false.neg.rate <- pl.true.v.called$fn / (pl.true.v.called$fn + pl.true.v.called$tp)
    pl.true.v.called$true.pos.rate <- pl.true.v.called$tp / (pl.true.v.called$tp + pl.true.v.called$fn)
    pl.true.v.called$true.neg.rate <- pl.true.v.called$tn / (pl.true.v.called$tn + pl.true.v.called$fp)
    pl.true.v.called$true.froh <- pl.true.v.called$tot.true.len/chrom.len
    write.csv(pl.true.v.called, paste0('bcftools_output/false_pos_neg_rates_',d,'_',h,'_PL.csv'),
              row.names = FALSE, quote = FALSE)


    ##### >>> 3B. BCFtools/ROH - GT (genotypes only) #####
    OUT <- NULL
    for(i in sort(unique(bcf.gt.res$id))){
      print(i)
      for(c in c(5,10,15,30,50)){
        ## get all called ROHs >= 100 kb
        called <- bcf.gt.res[bcf.gt.res$id == i & bcf.gt.res$covg == c,]
        tot.call.len <- sum(called[!duplicated(called$called.roh.id), 'length'])
        n.called <- length(unique(called$called.roh.id))

        ## get all true ROHs >= 100 kb
        true <- true.rohs[true.rohs$id == i & true.rohs$length >= 100000,]
        tot.true.len <- sum(true$length)
        n.true <- length(unique(true$true.roh.id))

        ## for each called ROH,
        TOT.TRUE <- 0
        for(r in unique(called$called.roh.id)){
          ## total up length that overlaps any >= 100 kb true ROH
          true.call.len <- sum(gt.overlap[gt.overlap$called.roh.id == r, 'true.len'])
          TOT.TRUE <- TOT.TRUE + true.call.len
        }
        save <- c(i, c, tot.call.len, n.called, tot.true.len, n.true, TOT.TRUE)
        OUT <- rbind(OUT, save)
      }
    }
    gt.true.v.called <- as.data.frame(OUT)
    colnames(gt.true.v.called) <- c('id','covg','total.call.len','n.called','tot.true.len','n.true','tp')
    for(c in 1:ncol(gt.true.v.called)){
      gt.true.v.called[,c] <- as.integer(gt.true.v.called[,c])
    }
    gt.true.v.called$fn <- gt.true.v.called$tot.true.len - gt.true.v.called$tp
    gt.true.v.called$fp <- gt.true.v.called$total.call.len - gt.true.v.called$tp
    gt.true.v.called$tn <- chrom.len - gt.true.v.called$tp - gt.true.v.called$fn - gt.true.v.called$fp
    gt.true.v.called$false.pos.rate <- gt.true.v.called$fp / (gt.true.v.called$fp + gt.true.v.called$tn)
    gt.true.v.called$false.neg.rate <- gt.true.v.called$fn / (gt.true.v.called$fn + gt.true.v.called$tp)
    gt.true.v.called$true.pos.rate <- gt.true.v.called$tp / (gt.true.v.called$tp + gt.true.v.called$fn)
    gt.true.v.called$true.neg.rate <- gt.true.v.called$tn / (gt.true.v.called$tn + gt.true.v.called$fp)
    gt.true.v.called$true.froh <- gt.true.v.called$tot.true.len/chrom.len
    write.csv(gt.true.v.called, paste0('bcftools_output/false_pos_neg_rates_',d,'_',h,'_GT.csv'),
              row.names = FALSE, quote = FALSE)


    ##### >>> 3C. PLINK #####
    OUT <- NULL
    for(i in unique(plink.res$id)){   ## ... for each individual
      print(i)
      for(c in unique(plink.res$covg)){
        ## get all called ROHs >= 100 kb
        called <- plink.res[plink.res$id == i & plink.res$covg == c,]
        tot.call.len <- sum(called[!duplicated(called$called.roh.id), 'length'])
        n.called <- length(unique(called$called.roh.id))

        ## get all true ROHs >= 100 kb
        true <- true.rohs[true.rohs$id == i & true.rohs$length >= 100000,]
        tot.true.len <- sum(true$length)   ## total length of all true ROHs > 100 kb in length
        n.true <- nrow(true)

        ## for each called ROH,
        TOT.TRUE <- 0 ## TOT.TRUE == tot.hit
        for(r in unique(called$called.roh.id)){
          ## total up length that overlaps any >= 100 kb true ROH
          true.call.len <- sum(plink.overlap[plink.overlap$called.roh.id == r, 'true.len'])
          TOT.TRUE <- TOT.TRUE + true.call.len
        }
        save <- c(unlist(called[1,c(1,5,6:13)], use.names = FALSE), tot.call.len, n.called, tot.true.len, n.true, TOT.TRUE)
        OUT <- rbind(OUT, save)
      }
    }
    plink.true.v.called <- as.data.frame(OUT)
    colnames(plink.true.v.called) <- c('id','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','total.call.len','n.called','tot.true.len','n.true','tp')
    for(c in 1:ncol(plink.true.v.called)){
      plink.true.v.called[,c] <- as.integer(plink.true.v.called[,c])
    }
    plink.true.v.called$fn <- plink.true.v.called$tot.true.len - plink.true.v.called$tp
    plink.true.v.called$fp <- plink.true.v.called$total.call.len - plink.true.v.called$tp
    plink.true.v.called$tn <- chrom.len - plink.true.v.called$tp - plink.true.v.called$fn - plink.true.v.called$fp
    plink.true.v.called$false.pos.rate <- plink.true.v.called$fp / (plink.true.v.called$fp + plink.true.v.called$tn)
    plink.true.v.called$false.neg.rate <- plink.true.v.called$fn / (plink.true.v.called$fn + plink.true.v.called$tp)
    plink.true.v.called$true.pos.rate <- plink.true.v.called$tp / (plink.true.v.called$tp + plink.true.v.called$fn)
    plink.true.v.called$true.neg.rate <- plink.true.v.called$tn / (plink.true.v.called$tn + plink.true.v.called$fp)
    plink.true.v.called$true.froh <- plink.true.v.called$tot.true.len/chrom.len
    write.csv(plink.true.v.called, paste0('plink_output/false_pos_neg_rates_',d,'_',h,'_PLINK.csv'),
              row.names = FALSE, quote = FALSE)

    ##### >>> 3D. Plotting false-pos and -neg rates for all 3 analyses (6 separate PDF files) #####
    ## all points plus mean +/- SE
    i <- 1
    while(i == 1){ ## a loop just to make all the separate plots write with 1 command
      text.size <- 1.25
      shrink <- 400 ## higher # here ==> narrower x-direction spread for points
      alph <- 0.2
      false.pos.y.lim <- 0.3

      ##### >>>>> False positives #####
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_PL_false_pos_rates.pdf'), width = 5, height = 6)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(pl.true.v.called$covg, pl.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage',
           ylab = 'False positive rate', col = 'transparent', main = 'Genotype likelihoods', xlim = c(0.5,5.5),
           cex.axis = text.size, cex.lab = text.size, xaxt = 'n')
        # axis(2, at = c(0, 0.1, 0.2, 0.3), labels = c('0.0','0.1','0.2','0.3'), cex.axis = text.size)
        axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
        x <- 1
        for(c in unique(pl.true.v.called$covg)){
          temp <- pl.true.v.called[pl.true.v.called$covg == c,]
          f <- sample(c(-100:100), length(unique(temp$id)))
          f <- f/shrink+x
          points(f, temp$false.pos.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
          points(x, mean(temp$false.pos.rate, na.rm = TRUE), pch = 16, col = pl.col)
          suppressWarnings(arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10),
                                  y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10),
                                  lwd = 2, col = pl.col, code=3, angle=90, length=0.1))
          x <- x+1
        }
      dev.off()

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_GT_false_pos_rates.pdf'), width = 5, height = 6)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(gt.true.v.called$covg, gt.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage',
           ylab = 'False positive rate', col = 'transparent', main = 'Genotypes only', xlim = c(0.5,5.5),
           cex.axis = text.size, cex.lab = text.size, xaxt = 'n')
      # axis(2, at = c(0, 0.1, 0.2, 0.3), labels = c('0.0','0.1','0.2','0.3'), cex.axis = text.size)
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      x <- 1
      for(c in unique(gt.true.v.called$covg)){
        temp <- gt.true.v.called[gt.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.pos.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.pos.rate, na.rm = TRUE), pch = 16, col = gt.col)
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10),
               y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10),
               lwd = 2, col = gt.col, code=3, angle=90, length=0.1)
        x <- x+1
      }
      dev.off()

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_PLINK_false_pos_rates.pdf'), width = 5, height = 6)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(plink.true.v.called$covg, plink.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage',
           ylab = 'False positive rate', col = 'transparent', main = 'PLINK', xlim = c(0.5,5.5),
           cex.axis = text.size, cex.lab = text.size, xaxt = 'n')
      # axis(2, at = c(0, 0.1, 0.2, 0.3), labels = c('0.0','0.1','0.2','0.3'), cex.axis = text.size)
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      x <- 1
      for(c in unique(plink.true.v.called$covg)){
        temp <- plink.true.v.called[plink.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.pos.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.pos.rate), pch = 16, col = plink.col)
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate) - sd(temp$false.pos.rate)/10), y1 = (mean(temp$false.pos.rate) + sd(temp$false.pos.rate)/10),
               lwd = 2, col = plink.col, code=3, angle=90, length=0.1)
        x <- x+1
      }
      dev.off()

      ##### >>>>> False negatives #####
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_PL_false_neg_rates.pdf'), width = 5, height = 6)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(pl.true.v.called$covg, pl.true.v.called$false.neg.rate, ylim = c(0,1), xlab = 'Coverage', ylab = 'False negative rate', col = 'transparent', main = 'Genotype likelihoods', xlim = c(0.5,5.5), cex.axis = text.size, cex.lab = text.size, xaxt = 'n')
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      x <- 1
      for(c in unique(pl.true.v.called$covg)){
        temp <- pl.true.v.called[pl.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.neg.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.neg.rate, na.rm = TRUE), pch = 16, col = pl.col)
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10),
               y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10),
               lwd = 2, col = pl.col, code=3, angle=90, length=0.1)
        x <- x+1
      }
      dev.off()

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_GT_false_neg_rates.pdf'), width = 5, height = 6)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(gt.true.v.called$covg, gt.true.v.called$false.neg.rate, ylim = c(0,1), xlab = 'Coverage', ylab = 'False negative rate', col = 'transparent', main = 'Genotypes only', xlim = c(0.5,5.5), cex.axis = text.size, cex.lab = text.size, xaxt = 'n')
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      x <- 1
      for(c in unique(gt.true.v.called$covg)){
        temp <- gt.true.v.called[gt.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.neg.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.neg.rate, na.rm = TRUE), pch = 16, col = gt.col)
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10),
               y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10),
               lwd = 2, col = gt.col, code=3, angle=90, length=0.1)
        x <- x+1
      }
      dev.off()

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_PLINK_false_neg_rates.pdf'), width = 5, height = 6)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(plink.true.v.called$covg, plink.true.v.called$false.neg.rate, ylim = c(0,1), xlab = 'Coverage', ylab = 'False negative rate', col = 'transparent', main = 'PLINK', xlim = c(0.5,5.5), xaxt = 'n', cex.axis = text.size, cex.lab = text.size)
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      x <- 1
      for(c in unique(plink.true.v.called$covg)){
        temp <- plink.true.v.called[plink.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.neg.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.neg.rate), pch = 16, col = plink.col)
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate) - sd(temp$false.neg.rate)/10), y1 = (mean(temp$false.neg.rate) + sd(temp$false.neg.rate)/10),
               lwd = 2, col = plink.col, code=3, angle=90, length=0.1)
        x <- x+1
      }
      dev.off()

      i <- i+1
    }

    ##### >>> 3E. True f(ROH) vs. false pos and neg rates per method #####
    ## unique colors for each method
    pl.pal <- colorRampPalette(c(lt.pl.col, max.covg.col))
    # pl.pal <- colorRampPalette(c(pl.col, max.covg.col))
    pl.cols <- pl.pal(5)
    col <- 1
    for(c in sort(unique(pl.true.v.called$covg))){
      pl.true.v.called[pl.true.v.called$covg == c, 'temp.colour'] <- pl.cols[col]
      col <- col+1
    }
    # gt.pal <- colorRampPalette(c(lt.gt.col, max.covg.col))
    gt.pal <- colorRampPalette(c(gt.col, max.covg.col))
    gt.cols <- gt.pal(5)
    col <- 1
    for(c in sort(unique(gt.true.v.called$covg))){
      gt.true.v.called[gt.true.v.called$covg == c, 'temp.colour'] <- gt.cols[col]
      col <- col+1
    }
    # plink.pal <- colorRampPalette(c(lt.plink.col, max.covg.col))
    plink.pal <- colorRampPalette(c(plink.col, max.covg.col))
    plink.cols <- plink.pal(5)
    col <- 1
    for(c in sort(unique(plink.true.v.called$covg))){
      plink.true.v.called[plink.true.v.called$covg == c, 'temp.colour'] <- plink.cols[col]
      col <- col+1
    }

    k <- 1
    while(k == 1){

      text.size <- 1.25
      alph <- 1
      poly.alph <- 0.4
      wid <- 2

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_true_fROH_vs_pos_neg_rates.pdf'), width = 7, height = 6)
      figure.ct <- figure.ct + 1
      par(mai = c(1.02,0.82,0.82,1.42))
      ## PL
      plot(pl.true.v.called$true.froh, pl.true.v.called$false.neg.rate, col = alpha(pl.true.v.called$temp.colour, alph), pch = 16,
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = 'False negative rate', main = 'Likelihoods', cex.lab = text.size, cex.axis = text.size)
      for(c in unique(pl.true.v.called$covg)){
        sub <- pl.true.v.called[pl.true.v.called$covg == c,]
        mod <- lm(sub$false.neg.rate ~ sub$true.froh)
        new.dat <- data.frame(true.froh=c(sort(sub$true.froh))) ## new data for prediction
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(-new.vals[,1]),]
        ## CI polygon
        # polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
        #               sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
        polygon(x = c(sort(sub$true.froh),
                      sort(sub$true.froh, decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = FALSE)), border = NA,
                col = alpha(pl.true.v.called[pl.true.v.called$covg == c, 'temp.colour'][1], poly.alph))
        # lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
        lines(sort(sub$true.froh), new.vals[,1], lwd = wid,
              col = pl.true.v.called[pl.true.v.called$covg == c, 'temp.colour'][1])
      }
      par(xpd = TRUE)
      legend('right', pch = 16, col = c(alpha(pl.cols[1], alph), alpha(pl.cols[2], alph), alpha(pl.cols[3], alph), alpha(pl.cols[4], alph), alpha(pl.cols[5], alph)), legend = c('5X','10X','15X','30X','50X'), inset = -0.2, cex = text.size)
      par(xpd = FALSE)

      plot(pl.true.v.called$true.froh, pl.true.v.called$false.pos.rate, col = alpha(pl.true.v.called$temp.colour, alph), pch = 16,
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = 'False positive rate', main = 'Likelihoods', cex.lab = text.size, cex.axis = text.size)
      for(c in unique(pl.true.v.called$covg)){
        sub <- pl.true.v.called[pl.true.v.called$covg == c,]
        mod <- lm(sub$false.pos.rate ~ sub$true.froh)
        new.dat <- data.frame(true.froh=c(sort(sub$true.froh))) ## new data for prediction
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(new.vals[,1]),]
        ## CI polygon
        # polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
        #               sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
        polygon(x = c(sort(sub$true.froh),
                      sort(sub$true.froh, decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
                col = alpha(pl.true.v.called[pl.true.v.called$covg == c, 'temp.colour'][1], poly.alph))
        # lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
        lines(sort(sub$true.froh), new.vals[,1], lwd = wid,
              col = pl.true.v.called[pl.true.v.called$covg == c, 'temp.colour'][1])
      }
      par(xpd = TRUE)
      legend('right', pch = 16, col = c(alpha(pl.cols[1], alph), alpha(pl.cols[2], alph), alpha(pl.cols[3], alph), alpha(pl.cols[4], alph), alpha(pl.cols[5], alph)), legend = c('5X','10X','15X','30X','50X'), inset = -0.2, cex = text.size)
      par(xpd = FALSE)

      ## GT
      plot(gt.true.v.called$true.froh, gt.true.v.called$false.neg.rate, col = alpha(gt.true.v.called$temp.colour, alph), pch = 16,
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = 'False negative rate', main = 'Genotypes', cex.lab = text.size, cex.axis = text.size)
      for(c in unique(gt.true.v.called$covg)){
        sub <- gt.true.v.called[gt.true.v.called$covg == c,]
        mod <- lm(sub$false.neg.rate ~ sub$true.froh)
        new.dat <- data.frame(true.froh=c(sort(sub$true.froh))) ## new data for prediction
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(-new.vals[,1]),]
        ## CI polygon
        # polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
        #               sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
        polygon(x = c(sort(sub$true.froh),
                      sort(sub$true.froh, decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = FALSE)), border = NA,
                col = alpha(gt.true.v.called[gt.true.v.called$covg == c, 'temp.colour'][1], poly.alph))
        # lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
        lines(sort(sub$true.froh), new.vals[,1], lwd = wid,
              col = gt.true.v.called[gt.true.v.called$covg == c, 'temp.colour'][1])
      }
      par(xpd = TRUE)
      legend('right', pch = 16, col = c(alpha(gt.cols[1], alph), alpha(gt.cols[2], alph), alpha(gt.cols[3], alph), alpha(gt.cols[4], alph), alpha(gt.cols[5], alph)), legend = c('5X','10X','15X','30X','50X'), inset = -0.2, cex = text.size)
      par(xpd = FALSE)

      plot(gt.true.v.called$true.froh, gt.true.v.called$false.pos.rate, col = alpha(gt.true.v.called$temp.colour, alph), pch = 16,
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = 'False positive rate', main = 'Genotypes', cex.lab = text.size, cex.axis = text.size)
      for(c in unique(gt.true.v.called$covg)){
        sub <- gt.true.v.called[gt.true.v.called$covg == c,]
        mod <- lm(sub$false.pos.rate ~ sub$true.froh)
        new.dat <- data.frame(true.froh=c(sort(sub$true.froh))) ## new data for prediction
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(new.vals[,1]),]
        ## CI polygon
        # polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
        #               sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
        polygon(x = c(sort(sub$true.froh),
                      sort(sub$true.froh, decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
                col = alpha(gt.true.v.called[gt.true.v.called$covg == c, 'temp.colour'][1], poly.alph))
        # lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
        lines(sort(sub$true.froh), new.vals[,1], lwd = wid,
              col = gt.true.v.called[gt.true.v.called$covg == c, 'temp.colour'][1])
      }
      par(xpd = TRUE)
      legend('right', pch = 16, col = c(alpha(gt.cols[1], alph), alpha(gt.cols[2], alph), alpha(gt.cols[3], alph), alpha(gt.cols[4], alph), alpha(gt.cols[5], alph)), legend = c('5X','10X','15X','30X','50X'), inset = -0.2, cex = text.size)
      par(xpd = FALSE)

      ## PLINK
      plot(plink.true.v.called$true.froh, plink.true.v.called$false.neg.rate, col = alpha(plink.true.v.called$temp.colour, alph), pch = 16,
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = 'False negative rate', main = 'Genotypes', cex.lab = text.size, cex.axis = text.size)
      for(c in unique(plink.true.v.called$covg)){
        sub <- plink.true.v.called[plink.true.v.called$covg == c,]
        mod <- lm(sub$false.neg.rate ~ sub$true.froh)
        new.dat <- data.frame(true.froh=c(sort(sub$true.froh))) ## new data for prediction
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(-new.vals[,1]),]
        ## CI polygon
        # polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
        #               sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
        polygon(x = c(sort(sub$true.froh),
                      sort(sub$true.froh, decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = FALSE)), border = NA,
                col = alpha(plink.true.v.called[plink.true.v.called$covg == c, 'temp.colour'][1], poly.alph))
        # lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
        lines(sort(sub$true.froh), new.vals[,1], lwd = wid,
              col = plink.true.v.called[plink.true.v.called$covg == c, 'temp.colour'][1])
      }
      par(xpd = TRUE)
      legend('right', pch = 16, col = c(alpha(plink.cols[1], alph), alpha(plink.cols[2], alph), alpha(plink.cols[3], alph), alpha(plink.cols[4], alph), alpha(plink.cols[5], alph)), legend = c('5X','10X','15X','30X','50X'), inset = -0.2, cex = text.size)
      par(xpd = FALSE)

      plot(plink.true.v.called$true.froh, plink.true.v.called$false.pos.rate, col = alpha(plink.true.v.called$temp.colour, alph), pch = 16,
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = 'False positive rate', main = 'PLINK', cex.lab = text.size, cex.axis = text.size)
      for(c in unique(plink.true.v.called$covg)){
        sub <- plink.true.v.called[plink.true.v.called$covg == c,]
        mod <- lm(sub$false.pos.rate ~ sub$true.froh)
        new.dat <- data.frame(true.froh=c(sort(sub$true.froh))) ## new data for prediction
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(new.vals[,1]),]
        polygon(x = c(sort(sub$true.froh),
                      sort(sub$true.froh, decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
                col = alpha(plink.true.v.called[plink.true.v.called$covg == c, 'temp.colour'][1], poly.alph))
        lines(sort(sub$true.froh), new.vals[,1], lwd = wid,
              col = plink.true.v.called[plink.true.v.called$covg == c, 'temp.colour'][1])
      }
      par(xpd = TRUE)
      legend('right', pch = 16, col = c(alpha(plink.cols[1], alph), alpha(plink.cols[2], alph), alpha(plink.cols[3], alph), alpha(plink.cols[4], alph), alpha(plink.cols[5], alph)), legend = c('5X','10X','15X','30X','50X'), inset = -0.2, cex = text.size)
      par(xpd = FALSE)

      k <- k+1
      dev.off()
    }


    ##### >>> 3F. Combining false - and false + into single plot for all 3 methods (2 plots) #####
    ## all points plus mean +/- SE
    i <- 1
    while(i == 1){ ## a loop just to make all the separate plots write with 1 command
      text.size <- 1.25
      shrink <- 2000 ## higher # here ==> narrower x-direction spread for points
      alph <- 0.2
      false.pos.y.lim <- 1
      diff <- 0.3 ## distance between GT and PL/PLINK
      OUT <- NULL ## for storing means

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_3_methods_false_pos_rates.pdf'), width = 8, height = 4.75)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      # plot(pl.true.v.called$covg, pl.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage', ylab = 'False positive rate', col = 'transparent', main = 'False positive rates', xlim = c(0.65,5.35), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
      plot(pl.true.v.called$covg, pl.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage', ylab = 'False positive rate', col = 'transparent', main = 'False positive rates', xlim = c(0.65,5.35), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
        axis(2, at = c(0, 0.1, 0.2, 0.3), labels = c('0.0','0.1','0.2','0.3'), cex.axis = text.size)
        axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
        x <- 1
        for(c in unique(pl.true.v.called$covg)){
          temp <- pl.true.v.called[pl.true.v.called$covg == c,]
          f <- sample(c(-100:100), length(unique(temp$id)))
          f <- f/shrink+x
          points(f, temp$false.pos.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
          points(x, mean(temp$false.pos.rate, na.rm = TRUE), pch = 16, col = pl.col)
          suppressWarnings(arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
                                  y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
                                  lwd = 2, col = pl.col, code=3, angle=90, length=0.1))
          OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 1))
          x <- x+1
        }

      x <- 1-diff
      for(c in unique(gt.true.v.called$covg)){
        temp <- gt.true.v.called[gt.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.pos.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.pos.rate, na.rm = TRUE), pch = 16, col = gt.col)
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
               y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
               lwd = 2, col = gt.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 2))
        x <- x+1
      }

      x <- 1+diff
      for(c in unique(plink.true.v.called$covg)){
        temp <- plink.true.v.called[plink.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.pos.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.pos.rate, na.rm = TRUE), pch = 16, col = plink.col)
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
               y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
               lwd = 2, col = plink.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 3))
        x <- x+1
      }

      ## turn on or off vertical dashed lines between coverage levels
      abline(v = 1.5, lty = 2, col = 'grey')
      abline(v = 2.5, lty = 2, col = 'grey')
      abline(v = 3.5, lty = 2, col = 'grey')
      abline(v = 4.5, lty = 2, col = 'grey')

      ## alternatively, connect within-method points with lines (this doesn't look good)
      # lines(OUT[c(1:5), 1], OUT[c(1:5), 2], col = pl.col)
      # lines(OUT[c(6:10), 1], OUT[c(6:10), 2], col = gt.col)
      # lines(OUT[c(11:15), 1], OUT[c(11:15), 2], col = plink.col)

      legend('topright', fill = c(gt.col, pl.col, plink.col), legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'),
             inset = 0.02, bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)

      dev.off()
      i <- i+1
    }

    i <- 1
    while(i == 1){ ## a loop just to make all the separate plots write with 1 command
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_3_methods_false_neg_rates.pdf'), width = 8, height = 4.75)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(pl.true.v.called$covg, pl.true.v.called$false.neg.rate, ylim = c(0, 1), xlab = 'Coverage', ylab = 'False negative rate', col = 'transparent', main = 'False negative rates', xlim = c(0.65,5.35), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = text.size)
      x <- 1
      for(c in unique(pl.true.v.called$covg)){
        temp <- pl.true.v.called[pl.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.neg.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.neg.rate, na.rm = TRUE), pch = 16, col = pl.col)

        ## SEs
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10),
        #        lwd = 2, col = pl.col, code=3, angle=90, length=0.1)

        ## SDs
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)),
        #        lwd = 2, col = pl.col, code=3, angle=90, length=0.1)

        ## 95% CIs
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - (1.96*sd(temp$false.neg.rate, na.rm = TRUE)/10)),
               y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + (1.96*sd(temp$false.neg.rate, na.rm = TRUE)/10)),
               lwd = 2, col = pl.col, code=3, angle=90, length=0.1)

        OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 1))
        x <- x+1

      }

      x <- 1-diff
      for(c in unique(gt.true.v.called$covg)){
        temp <- gt.true.v.called[gt.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.neg.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.neg.rate, na.rm = TRUE), pch = 16, col = gt.col)

        ## SEs
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10),
        #        lwd = 2, col = gt.col, code=3, angle=90, length=0.1)

        ## SDs
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)),
        #        lwd = 2, col = gt.col, code=3, angle=90, length=0.1)

        ## 95% CIs
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - (1.96*sd(temp$false.neg.rate, na.rm = TRUE)/10)),
               y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + (1.96*sd(temp$false.neg.rate, na.rm = TRUE)/10)),
               lwd = 2, col = gt.col, code=3, angle=90, length=0.1)

        OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 2))
        x <- x+1
      }

      x <- 1+diff
      for(c in unique(plink.true.v.called$covg)){
        temp <- plink.true.v.called[plink.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        points(f, temp$false.neg.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
        points(x, mean(temp$false.neg.rate, na.rm = TRUE), pch = 16, col = plink.col)

        ## SEs
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10),
        #        lwd = 2, col = plink.col, code=3, angle=90, length=0.1)

        ## SDs
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)),
        #        lwd = 2, col = plink.col, code=3, angle=90, length=0.1)

        ## 95% CIs
        arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - (1.96*sd(temp$false.neg.rate, na.rm = TRUE)/10)),
               y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + (1.96*sd(temp$false.neg.rate, na.rm = TRUE)/10)),
               lwd = 2, col = plink.col, code=3, angle=90, length=0.1)

        OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 3))
        x <- x+1
      }

      ## turn on or off vertical dashed lines between coverage levels
      abline(v = 1.5, lty = 2, col = 'grey')
      abline(v = 2.5, lty = 2, col = 'grey')
      abline(v = 3.5, lty = 2, col = 'grey')
      abline(v = 4.5, lty = 2, col = 'grey')

      ## alternatively, connect within-method points with lines (this doesn't look good)
      # lines(OUT[c(1:5), 1], OUT[c(1:5), 2], col = pl.col)
      # lines(OUT[c(6:10), 1], OUT[c(6:10), 2], col = gt.col)
      # lines(OUT[c(11:15), 1], OUT[c(11:15), 2], col = plink.col)

      # legend('topright', fill = c(gt.col, pl.col, plink.col), legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'),
      #        inset = 0.02, bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)

      dev.off()
      i <- i+1
    }

    ##### >>> 3FII. Combining false - and false + into single plot for all 3 methods (2 custom boxplots) #####
    i <- 1
    while(i == 1){ ## a loop just to make all the separate plots write with 1 command
      text.size <- 1.25
      shrink <- 2000 ## higher # here ==> narrower x-direction spread for points
      alph <- 0.2
      poly.alph <- 0.4
      false.pos.y.lim <- 0.2
      percentile <- 50
      low.q <- (100-percentile)/100/2
      hi.q <- 1 - (100-percentile)/100/2
      diff <- 0.3 ## distance between GT and PL/PLINK
      OUT <- NULL ## for storing means

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_3_methods_false_pos_rates.pdf'), width = 8, height = 4.75)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      # plot(pl.true.v.called$covg, pl.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage', ylab = 'False positive rate', col = 'transparent', main = 'False positive rates', xlim = c(0.65,5.35), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
      plot(pl.true.v.called$covg, pl.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage', ylab = 'False positive rate', col = 'transparent', main = paste0('False positive rates - ',d), xlim = c(0.65,5.35), 
           cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
      axis(2, at = c(0, 0.1, 0.2, 0.3), labels = c('0.0','0.1','0.2','0.3'), cex.axis = text.size)
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      x <- 1
      for(c in unique(pl.true.v.called$covg)){
        temp <- pl.true.v.called[pl.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x

        polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                     quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                     quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                     quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2]),
                col = alpha(pl.col, poly.alph), border = NA)

        lines(x = c(x-diff/2, x+diff/2), y = c(median(temp$false.pos.rate, na.rm = TRUE), median(temp$false.pos.rate, na.rm = TRUE))
              , pch = 16, col = pl.col, lwd = 2)
        points(f, temp$false.pos.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
        #        y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
        #        lwd = 2, col = pl.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 1))
        x <- x+1

      }

      x <- 1-diff
      for(c in unique(gt.true.v.called$covg)){
        temp <- gt.true.v.called[gt.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x

        polygon(x = c(x-diff/2.2, x+diff/2.2, x+diff/2.2, x-diff/2.2), y = c(quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                             quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2]),
                col = alpha(gt.cols[3], poly.alph), border = NA)
        lines(x = c(x-diff/2.2, x+diff/2.2), y = c(median(temp$false.pos.rate, na.rm = TRUE), median(temp$false.pos.rate, na.rm = TRUE))
              , pch = 16, col = gt.col, lwd = 2)
        points(f, temp$false.pos.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
        #        y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
        #        lwd = 2, col = gt.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 1))
        x <- x+1
      }

      x <- 1+diff
      for(c in unique(plink.true.v.called$covg)){
        temp <- plink.true.v.called[plink.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        polygon(x = c(x-diff/2.2, x+diff/2.2, x+diff/2.2, x-diff/2.2), y = c(quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                             quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2]),
                col = alpha(plink.col, poly.alph), border = NA)
        lines(x = c(x-diff/2.2, x+diff/2.2), y = c(median(temp$false.pos.rate, na.rm = TRUE), median(temp$false.pos.rate, na.rm = TRUE))
              , pch = 16, col = plink.col, lwd = 2)
        points(f, temp$false.pos.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.pos.rate, na.rm = TRUE) - sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
        #        y1 = (mean(temp$false.pos.rate, na.rm = TRUE) + sd(temp$false.pos.rate, na.rm = TRUE)/10*1.96),
        #        lwd = 2, col = plink.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 1))
        x <- x+1
      }

      ## turn on or off vertical dashed lines between coverage levels
      abline(v = 1.5, lty = 2, col = 'grey')
      abline(v = 2.5, lty = 2, col = 'grey')
      abline(v = 3.5, lty = 2, col = 'grey')
      abline(v = 4.5, lty = 2, col = 'grey')

      ## alternatively, connect within-method points with lines (this doesn't look good)
      # lines(OUT[c(1:5), 1], OUT[c(1:5), 2], col = pl.col)
      # lines(OUT[c(6:10), 1], OUT[c(6:10), 2], col = gt.col)
      # lines(OUT[c(11:15), 1], OUT[c(11:15), 2], col = plink.col)

      legend('topright', fill = c(gt.col, pl.col, plink.col), legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'),
             inset = 0.02, bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)

      dev.off()

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_3_methods_false_neg_rates.pdf'), width = 8, height = 4.75)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1,6.1,4.1,2.1))
      plot(pl.true.v.called$covg, pl.true.v.called$false.neg.rate, ylim = c(0, 1), xlab = 'Coverage', ylab = 'False negative rate', col = 'transparent', main = paste0('False negative rates - ',d), xlim = c(0.65,5.35), 
           cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
      axis(1, at = c(1:5), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = text.size)
      x <- 1
      for(c in unique(pl.true.v.called$covg)){
        temp <- pl.true.v.called[pl.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x

        polygon(x = c(x-diff/2.2, x+diff/2.2, x+diff/2.2, x-diff/2.2), y = c(quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2]),
                col = alpha(pl.cols[3], poly.alph), border = NA)

        lines(x = c(x-diff/2.2, x+diff/2.2), y = c(median(temp$false.neg.rate, na.rm = TRUE), median(temp$false.neg.rate, na.rm = TRUE))
              , pch = 16, col = pl.col, lwd = 2)
        points(f, temp$false.neg.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10*1.96),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10*1.96),
        #        lwd = 2, col = pl.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 1))
        x <- x+1

      }

      x <- 1-diff
      for(c in unique(gt.true.v.called$covg)){
        temp <- gt.true.v.called[gt.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x

        polygon(x = c(x-diff/2.2, x+diff/2.2, x+diff/2.2, x-diff/2.2), y = c(quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2]),
                col = alpha(gt.cols[3], poly.alph), border = NA)
        lines(x = c(x-diff/2.2, x+diff/2.2), y = c(median(temp$false.neg.rate, na.rm = TRUE), median(temp$false.neg.rate, na.rm = TRUE))
              , pch = 16, col = gt.col, lwd = 2)
        points(f, temp$false.neg.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10*1.96),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10*1.96),
        #        lwd = 2, col = gt.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 1))
        x <- x+1
      }

      x <- 1+diff
      for(c in unique(plink.true.v.called$covg)){
        temp <- plink.true.v.called[plink.true.v.called$covg == c,]
        f <- sample(c(-100:100), length(unique(temp$id)))
        f <- f/shrink+x
        polygon(x = c(x-diff/2.2, x+diff/2.2, x+diff/2.2, x-diff/2.2), y = c(quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                             quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2]),
                col = alpha(plink.col, poly.alph), border = NA)
        lines(x = c(x-diff/2.2, x+diff/2.2), y = c(median(temp$false.neg.rate, na.rm = TRUE), median(temp$false.neg.rate, na.rm = TRUE))
              , pch = 16, col = plink.col, lwd = 2)
        points(f, temp$false.neg.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
        # arrows(x0 = x, x1 = x, y0 = (mean(temp$false.neg.rate, na.rm = TRUE) - sd(temp$false.neg.rate, na.rm = TRUE)/10*1.96),
        #        y1 = (mean(temp$false.neg.rate, na.rm = TRUE) + sd(temp$false.neg.rate, na.rm = TRUE)/10*1.96),
        #        lwd = 2, col = plink.col, code=3, angle=90, length=0.1)
        OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 1))
        x <- x+1
      }

      ## turn on or off vertical dashed lines between coverage levels
      abline(v = 1.5, lty = 2, col = 'grey')
      abline(v = 2.5, lty = 2, col = 'grey')
      abline(v = 3.5, lty = 2, col = 'grey')
      abline(v = 4.5, lty = 2, col = 'grey')

      ## alternatively, connect within-method points with lines (this doesn't look good)
      # lines(OUT[c(1:5), 1], OUT[c(1:5), 2], col = pl.col)
      # lines(OUT[c(6:10), 1], OUT[c(6:10), 2], col = gt.col)
      # lines(OUT[c(11:15), 1], OUT[c(11:15), 2], col = plink.col)

      # legend('topright', fill = c(gt.col, pl.col, plink.col), legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'),
      #        inset = 0.02, bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)

      dev.off()
      i <- i+1
    }

    ##### >>> 3G. Organizing false neg/pos rates into table #####
    OUT <- NULL
    for(c in sort(unique(plink.true.v.called$covg))){
      gt.median.pos <- median(gt.true.v.called[gt.true.v.called$covg == c, 'false.pos.rate'])
      gt.mean.pos <- mean(gt.true.v.called[gt.true.v.called$covg == c, 'false.pos.rate'])
      gt.sd.pos <- sd(gt.true.v.called[gt.true.v.called$covg == c, 'false.pos.rate'])
      gt.cv.pos <- gt.sd.pos/gt.mean.pos*100
      gt.median.neg <- median(gt.true.v.called[gt.true.v.called$covg == c, 'false.neg.rate'])
      gt.mean.neg <- mean(gt.true.v.called[gt.true.v.called$covg == c, 'false.neg.rate'])
      gt.sd.neg <- sd(gt.true.v.called[gt.true.v.called$covg == c, 'false.neg.rate'])
      gt.cv.neg <- gt.sd.neg/gt.mean.neg*100

      pl.median.pos <- median(pl.true.v.called[pl.true.v.called$covg == c, 'false.pos.rate'])
      pl.mean.pos <- mean(pl.true.v.called[pl.true.v.called$covg == c, 'false.pos.rate'])
      pl.sd.pos <- sd(pl.true.v.called[pl.true.v.called$covg == c, 'false.pos.rate'])
      pl.cv.pos <- pl.sd.pos/pl.mean.pos*100
      pl.median.neg <- median(pl.true.v.called[pl.true.v.called$covg == c, 'false.neg.rate'])
      pl.mean.neg <- mean(pl.true.v.called[pl.true.v.called$covg == c, 'false.neg.rate'])
      pl.sd.neg <- sd(pl.true.v.called[pl.true.v.called$covg == c, 'false.neg.rate'])
      pl.cv.neg <- pl.sd.neg/pl.mean.neg*100

      plink.median.pos <- median(plink.true.v.called[plink.true.v.called$covg == c, 'false.pos.rate'])
      plink.mean.pos <- mean(plink.true.v.called[plink.true.v.called$covg == c, 'false.pos.rate'])
      plink.sd.pos <- sd(plink.true.v.called[plink.true.v.called$covg == c, 'false.pos.rate'])
      plink.cv.pos <- plink.sd.pos/plink.mean.pos*100
      plink.median.neg <- median(plink.true.v.called[plink.true.v.called$covg == c, 'false.neg.rate'])
      plink.mean.neg <- mean(plink.true.v.called[plink.true.v.called$covg == c, 'false.neg.rate'])
      plink.sd.neg <- sd(plink.true.v.called[plink.true.v.called$covg == c, 'false.neg.rate'])
      plink.cv.neg <- plink.sd.neg/plink.mean.neg*100

      save <- rbind(c(c, 'Genotypes', gt.median.pos, gt.mean.pos, gt.sd.pos, gt.cv.pos, gt.median.neg, gt.mean.neg, gt.sd.neg, gt.cv.neg),
                    c(c, 'Likelihoods', pl.median.pos, pl.mean.pos, pl.sd.pos, pl.cv.pos, pl.median.neg, pl.mean.neg, pl.sd.neg, pl.cv.neg),
                    c(c, 'PLINK', plink.median.pos, plink.mean.pos, plink.sd.pos, plink.cv.pos, plink.median.neg, plink.mean.neg, plink.sd.neg, plink.cv.neg))
      OUT <- rbind(OUT, save)
    }
    colnames(OUT) <- c('covg','method','median.fp','mean.fp','sd.fp','cv.fp','median.fn','mean.fn','sd.fn','cv.fn')
    write.csv(OUT, paste0('../data/3_methods_results/',d,'_false_neg_pos_rate_stats.csv'), row.names = FALSE)

    ##### 4. Plotting true vs. called f(ROH) values #####
    ## assign colors based on coverage level
    ## same colors across methods
    # pal <- colorRampPalette(c(min.covg.col, max.covg.col))
    # cols <- pal(5)
    # col <- 1
    # for(c in sort(unique(froh.stats$covg))){
    #   froh.stats[froh.stats$covg == c, 'pl.temp.col'] <- cols[col]
    #   froh.stats[froh.stats$covg == c, 'gt.temp.col'] <- cols[col]
    #   froh.stats[froh.stats$covg == c, 'plink.temp.col'] <- cols[col]
    #   col <- col+1
    # }

    ## unique colors for each method
    k <- 1
    while(k == 1){

      pl.pal <- colorRampPalette(c(lt.pl.col, max.covg.col))
      # pl.pal <- colorRampPalette(c(pl.col, max.covg.col))
      pl.cols <- pl.pal(5)
      col <- 1
      for(c in sort(unique(froh.stats$covg))){
        froh.stats[froh.stats$covg == c, 'pl.temp.col'] <- pl.cols[col]
        col <- col+1
      }
      # gt.pal <- colorRampPalette(c(lt.gt.col, max.covg.col))
      gt.pal <- colorRampPalette(c(gt.col, max.covg.col))
      gt.cols <- gt.pal(5)
      col <- 1
      for(c in sort(unique(froh.stats$covg))){
        froh.stats[froh.stats$covg == c, 'gt.temp.col'] <- gt.cols[col]
        col <- col+1
      }
      # plink.pal <- colorRampPalette(c(lt.plink.col, max.covg.col))
      plink.pal <- colorRampPalette(c(plink.col, max.covg.col))
      plink.cols <- plink.pal(5)
      col <- 1
      for(c in sort(unique(froh.stats$covg))){
        froh.stats[froh.stats$covg == c, 'plink.temp.col'] <- plink.cols[col]
        col <- col+1
      }

      wid <- 2
      alph <- 0.4
      txt.size <- 1.5
      poly.alph = 0.25

      new.dat <- data.frame(true.froh=c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]))) ## new data for prediction

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_true_vs_called_fROH.pdf'), width = 4, height = 12)
      par(mar = c(5.1, 4.6, 4.1, 2.1))
      figure.ct <- figure.ct + 1
      par(mfrow = c(3,1))
      
      ## calculate y-lim values based on 95% CIs
      lo.y <- min(froh.stats$pl.froh, froh.stats$gt.froh, froh.stats$plink.froh)
      hi.y <- max(froh.stats$pl.froh, froh.stats$gt.froh, froh.stats$plink.froh)
      for(c in unique(froh.stats$covg)){
        mod <- lm(froh.stats[froh.stats$covg == c, 'pl.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence', level = 0.95)
        new.vals <- new.vals[order(new.vals[,1]),]
        lo.y <- min(lo.y, min(new.vals[,2]))
        hi.y <- max(hi.y, max(new.vals[,3]))
      }
      for(c in unique(froh.stats$covg)){
        mod <- lm(froh.stats[froh.stats$covg == c, 'gt.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence', level = 0.95)
        new.vals <- new.vals[order(new.vals[,1]),]
        lo.y <- min(lo.y, min(new.vals[,2]))
        hi.y <- max(hi.y, max(new.vals[,3]))
      }
      for(c in unique(froh.stats$covg)){
        mod <- lm(froh.stats[froh.stats$covg == c, 'plink.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence', level = 0.95)
        new.vals <- new.vals[order(new.vals[,1]),]
        lo.y <- min(lo.y, min(new.vals[,2]))
        hi.y <- max(hi.y, max(new.vals[,3]))
      }
      
      ## GT
      plot(froh.stats$true.froh, froh.stats$gt.froh, pch = 16, col = 'transparent', main = 'Genotypes only',
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = substitute(paste('Called ',italic('F')[ROH])), 
           cex.axis = txt.size, cex.lab = txt.size,
           # xlim = c(0,1), ylim = c(0, 1))
           xlim = c(min(froh.stats$true.froh), max(froh.stats$true.froh)),
           ylim = c(lo.y, hi.y))
      abline(0,1, lty = 2)
      points(froh.stats$true.froh, froh.stats$gt.froh, pch = 16, col = alpha(froh.stats$gt.temp.col, alph))
      for(c in unique(froh.stats$covg)){
        mod <- lm(froh.stats[froh.stats$covg == c, 'gt.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(new.vals[,1]),]
        ## CI polygon
        polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
                      sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
                col = alpha(froh.stats[froh.stats$covg == c, 'gt.temp.col'][1], poly.alph))
        lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
              col = froh.stats[froh.stats$covg == c, 'gt.temp.col'][1])
      }
      if(d == 'large-1000'){
        legend('topleft', col = gt.cols, legend = c('5X','10X','15X','30X','50X'), lwd = wid, inset = 0, bty = 'n', cex = txt.size)
      }
      
      ## PL
      plot(froh.stats$true.froh, froh.stats$pl.froh, pch = 16, col = 'transparent', main = 'Genotype likelihoods',
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = substitute(paste('Called ',italic('F')[ROH])), 
           cex.axis = txt.size, cex.lab = txt.size,
           # xlim = c(0,1), ylim = c(0, 1))
           xlim = c(min(froh.stats$true.froh), max(froh.stats$true.froh)),
           ylim = c(lo.y, hi.y))
      abline(0,1, lty = 2)
      points(froh.stats$true.froh, froh.stats$pl.froh, pch = 16, col = alpha(froh.stats$pl.temp.col, alph))
      for(c in unique(froh.stats$covg)){
        mod <- lm(froh.stats[froh.stats$covg == c, 'pl.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence', level = 0.95)
        new.vals <- new.vals[order(new.vals[,1]),]

        ## CI polygon
        polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
                      sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
                col = alpha(froh.stats[froh.stats$covg == c, 'pl.temp.col'][1], poly.alph))
        lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
              col = froh.stats[froh.stats$covg == c, 'pl.temp.col'][1])

      }
      if(d == 'large-1000'){
        legend('topleft', col = pl.cols, legend = c('5X','10X','15X','30X','50X'), lwd = wid, inset = 0, bty = 'n', cex = txt.size)
      }
      
      ## PLINK
      plot(froh.stats$true.froh, froh.stats$plink.froh, pch = 16, col = 'transparent', main = 'PLINK',
           xlab = substitute(paste('True ',italic('F')[ROH])), ylab = substitute(paste('Called ',italic('F')[ROH])), 
           cex.axis = txt.size, cex.lab = txt.size,
           # xlim = c(0,1), ylim = c(0, 1))
           xlim = c(min(froh.stats$true.froh), max(froh.stats$true.froh)),
           ylim = c(lo.y, hi.y))
      abline(0,1, lty = 2)
      points(froh.stats$true.froh, froh.stats$plink.froh, pch = 16, col = alpha(froh.stats$plink.temp.col, alph))
      for(c in unique(froh.stats$covg)){
        mod <- lm(froh.stats[froh.stats$covg == c, 'plink.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
        new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
        new.vals <- new.vals[order(new.vals[,1]),]
        ## CI polygon
        polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
                      sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
                y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
                col = alpha(froh.stats[froh.stats$covg == c, 'plink.temp.col'][1], poly.alph))
        lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
              col = froh.stats[froh.stats$covg == c, 'plink.temp.col'][1])
      }
      if(d == 'large-1000'){
        legend('topleft', col = plink.cols, legend = c('5X','10X','15X','30X','50X'), lwd = wid, inset = 0, bty = 'n', cex = txt.size)
      }
      dev.off()
      k <- k+1
    }

    new.dat <- data.frame(true.froh=c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]))) ## new data for prediction
    ## all on single plot for stats viz
    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_all_methods_true_vs_called_fROH.pdf'), width = 7, height = 7)
    figure.ct <- figure.ct + 1
    plot(c(froh.stats$true.froh, froh.stats$true.froh, froh.stats$true.froh),
         c(froh.stats$pl.froh, froh.stats$gt.froh, froh.stats$plink.froh), pch = 16, col = 'transparent', main = 'All methods',
         xlab = substitute(paste('True ',italic('F')[ROH])), ylab = substitute(paste('Called ',italic('F')[ROH])), cex.axis = txt.size, cex.lab = txt.size)
    abline(0,1, lty = 2)
    points(froh.stats$true.froh, froh.stats$pl.froh, pch = 16, col = alpha(froh.stats$pl.temp.col, alph))
    for(c in unique(froh.stats$covg)){
      mod <- lm(froh.stats[froh.stats$covg == c, 'pl.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
      # suppressWarnings(clipplot(abline(mod, col = froh.stats[froh.stats$covg == c, 'pl.temp.col'][1]),
      #                           xlim = c(min(froh.stats[froh.stats$covg == c, 'true.froh']),
      #                                    max(froh.stats[froh.stats$covg == c, 'true.froh']))))
      new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
      new.vals <- new.vals[order(new.vals[,1]),]
      ## CI polygon
      polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
                    sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
              y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
              col = alpha(froh.stats[froh.stats$covg == c, 'pl.temp.col'][1], poly.alph))
      lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
            col = froh.stats[froh.stats$covg == c, 'pl.temp.col'][1])
    }
    points(froh.stats$true.froh, froh.stats$gt.froh, pch = 16, col = alpha(froh.stats$gt.temp.col, alph))
    for(c in unique(froh.stats$covg)){
      mod <- lm(froh.stats[froh.stats$covg == c, 'gt.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
      # suppressWarnings(clipplot(abline(mod, col = froh.stats[froh.stats$covg == c, 'pl.temp.col'][1]),
      #                           xlim = c(min(froh.stats[froh.stats$covg == c, 'true.froh']),
      #                                    max(froh.stats[froh.stats$covg == c, 'true.froh']))))
      new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
      new.vals <- new.vals[order(new.vals[,1]),]
      ## CI polygon
      polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
                    sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
              y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
              col = alpha(froh.stats[froh.stats$covg == c, 'gt.temp.col'][1], poly.alph))
      lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
            col = froh.stats[froh.stats$covg == c, 'gt.temp.col'][1])
    }
    points(froh.stats$true.froh, froh.stats$plink.froh, pch = 16, col = alpha(froh.stats$plink.temp.col, alph))
    for(c in unique(froh.stats$covg)){
      mod <- lm(froh.stats[froh.stats$covg == c, 'plink.froh'] ~ froh.stats[froh.stats$covg == c, 'true.froh'])
      # suppressWarnings(clipplot(abline(mod, col = froh.stats[froh.stats$covg == c, 'pl.temp.col'][1]),
      #                           xlim = c(min(froh.stats[froh.stats$covg == c, 'true.froh']),
      #                                    max(froh.stats[froh.stats$covg == c, 'true.froh']))))
      new.vals <- predict(mod, newdata = new.dat, interval = 'confidence')
      new.vals <- new.vals[order(new.vals[,1]),]
      ## CI polygon
      polygon(x = c(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]),
                    sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3], decreasing = TRUE)),
              y = c(new.vals[,2], sort(new.vals[,3], decreasing = TRUE)), border = NA,
              col = alpha(froh.stats[froh.stats$covg == c, 'plink.temp.col'][1], poly.alph))
      lines(sort(froh.stats[!duplicated(froh.stats[,c(1,3)]), 3]), new.vals[,1], lwd = wid,
            col = froh.stats[froh.stats$covg == c, 'plink.temp.col'][1])
    }
    legend('topleft', col = c(pl.col, gt.col, plink.col), lwd = wid, legend = c('PL','GT','PLINK'),
           bty = 'n', inset = 0.02, cex = txt.size)
    dev.off()

    ##### >>> 4A. Statistics: comparing true vs. called f(ROH) slopes and intercepts #####
    ## reformat data for lm
    dat <- rbind(cbind(froh.stats[, c('id','covg','true.froh')], call.froh = froh.stats$pl.froh,
                       method = rep('PL', length(unique(froh.stats$id)))),
                 cbind(froh.stats[, c('id','covg','true.froh')], call.froh = froh.stats$gt.froh,
                       method = rep('GT', length(unique(froh.stats$id)))),
                 cbind(froh.stats[, c('id','covg','true.froh')], call.froh = froh.stats$plink.froh,
                       method = rep('PLINK', length(unique(froh.stats$id)))))
    for(c in 1:4){
      dat[,c] <- as.numeric(dat[,c])
    }

    ## check assumptions and some stats
    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_truefROH_vs_callfROH_lm_assumptions.pdf'), width = 10, height = 10)
    figure.ct <- figure.ct + 1
    par(mfrow = c(2,2))
    OUT <- NULL ## for saving statsy info
    for(m in unique(dat$method)){
      for(c in sort(unique(dat$covg))){
        print(paste0(m,' - ',c))
        sub <- dat[dat$method == m & dat$covg == c,]
        mod <- lm(sub$call.froh ~ sub$true.froh)
        # mod <- lm(log(sub$call.froh + 1) ~ sub$true.froh) ## trying log x-form to deal with heteroscedasticity in PL models (didn't help)
        int.lo.lim <- confint(mod, level = 0.95)[1,1]
        int.up.lim <- confint(mod, level = 0.95)[1,2]
        slope.lo.lim <- confint(mod, level = 0.95)[2,1]
        slope.up.lim <- confint(mod, level = 0.95)[2,2]
        int <- mod$coefficients[1]
        slope <- mod$coefficients[2]
        save <- c(m, c, summary(mod)$adj.r.squared, summary(mod)$coefficients[2,4], int.lo.lim, int, int.up.lim, slope.lo.lim, slope, slope.up.lim)
        OUT <- rbind(OUT, save)

        ## test for normality of residuals (none print, so all good here)
        # if(shapiro.test(residuals(mod))[2] < 0.05){
        #   print(paste0('Residuals not normal: ',m,' - ',c,'X'))
        # }
        plot(mod, main = paste0(m,' - ',c,'X'))

        ## test whether intercept differs from 0 using 95% CIs
        if(confint(mod, level = 0.95)[1,1] <= 0 & confint(mod, level = 0.95)[1,2] >= 0){
          print(paste0(m,' - ',c,'X - Intercept == 0'))
        }
        if((confint(mod, level = 0.95)[1,1] <= 0 & confint(mod, level = 0.95)[1,2] <= 0) |
           (confint(mod, level = 0.95)[1,1] >= 0 & confint(mod, level = 0.95)[1,2] >= 0)){
          print(paste0(m,' - ',c,'X - Intercept != 0'))
        }

        ## test whether slopes differ from 1 using 95% CIs
        if(confint(mod, level = 0.95)[2,1] <= 1 & confint(mod, level = 0.95)[2,2] >= 1){
          print(paste0(m,' - ',c,'X - Slope == 1'))
        }
        if((confint(mod, level = 0.95)[2,1] <= 1 & confint(mod, level = 0.95)[2,2] <= 1) |
           (confint(mod, level = 0.95)[2,1] >= 1 & confint(mod, level = 0.95)[2,2] >= 1)){
          print(paste0(m,' - ',c,'X - Slope != 1'))
        }
        confint(mod, level = 0.95)
        # print(summary(mod))
      }
    }
    dev.off()
    colnames(OUT) <- c('method','covg','adj.r2','p.val','int.lo.lim','int','int.up.lim','slope.lo.lim','slope','slope.up.lim')
    write.csv(OUT, paste0('../data/3_methods_results/',d,'_true_v_callfROH_stats.csv'), row.names = FALSE)

    ##### >>> 4B. Plotting true f(ROH) vs. called - true f(ROH) #####
    dat$diff <- dat$call.froh - dat$true.froh
    alph <- 0.5
    y.min <- -0.7
    y.max <- 0.7

    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_true_vs_diff.pdf'), width = 6, height = 5)
    figure.ct <- figure.ct + 1

    ## PL
    sub <- dat[dat$method == 'PL',]
    plot(sub$true.froh, sub$diff, col = 'transparent', xlab = substitute(paste('True ',italic('F')[ROH])), 
         ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])),
         main = 'Likelihoods', ylim = c(y.min, y.max))
    abline(h = 0, lty = 2, col = 'grey')
    x <- 1
    for(c in sort(unique(sub$covg))){
      temp <- sub[sub$covg == c,]
      points(temp$true.froh, temp$diff, col = alpha(pl.cols[x], alph), pch = 16)
      abline(lm(temp$dif ~ temp$true.froh), col = pl.cols[x])
      x <- x+1
    }
    legend('topleft', col = alpha(pl.cols, alph), legend = c('5X','10X','15X','30X','50X'), pch = 16, bg = 'white')

    ## GT
    sub <- dat[dat$method == 'GT',]
    plot(sub$true.froh, sub$diff, col = 'transparent', xlab = substitute(paste('True ',italic('F')[ROH])), 
         ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])),
         main = 'Genotypes', ylim = c(y.min, y.max))
    abline(h = 0, lty = 2, col = 'grey')
    x <- 1
    for(c in sort(unique(sub$covg))){
      temp <- sub[sub$covg == c,]
      points(temp$true.froh, temp$diff, col = alpha(gt.cols[x], alph), pch = 16)
      abline(lm(temp$dif ~ temp$true.froh), col = gt.cols[x])
      x <- x+1
    }
    legend('topright', col = alpha(gt.cols, alph), legend = c('5X','10X','15X','30X','50X'), pch = 16, bg = 'white')

    ## PLINK
    sub <- dat[dat$method == 'PLINK',]
    plot(sub$true.froh, sub$diff, col = 'transparent', xlab = substitute(paste('True ',italic('F')[ROH])), 
         ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])),
         main = 'PLINK', ylim = c(y.min, y.max))
    abline(h = 0, lty = 2, col = 'grey')
    x <- 1
    for(c in sort(unique(sub$covg))){
      temp <- sub[sub$covg == c,]
      points(temp$true.froh, temp$diff, col = alpha(plink.cols[x], alph), pch = 16)
      abline(lm(temp$dif ~ temp$true.froh), col = plink.cols[x])
      x <- x+1
    }
    legend('topright', col = alpha(plink.cols, alph), legend = c('5X','10X','15X','30X','50X'), pch = 16, bg = 'white')
    dev.off()


    ##### 5. Plotting f(ROH) across coverage levels, individual lines ##### line 3 point 5
    alph <- 0.3  ## 1 for n = 30, 0.5 for n = 100
    bg.alph <- 0.2 ## for plot with mean lines +/- CIs (was previously 0.10)
    pt.alph <- 0.3
    lwd <- 1.5
    pt.size <- 1.5
    err.wid <- 0.05
    true.col <- 'black'

    txt.size <- 1.25
    ymin <- 0
    ymax <- 1
    xmin <- 1.85
    xmax <- 6.15
    n <- length(unique(froh.stats$id))

    OUT <- NULL
    for(c in sort(unique(froh.stats$covg))){
      save <- c(c,
                mean(froh.stats[froh.stats$covg == c, 'true.froh']),
                sd(froh.stats[froh.stats$covg == c, 'true.froh'])/sqrt(n),
                mean(froh.stats[froh.stats$covg == c, 'pl.froh']),
                sd(froh.stats[froh.stats$covg == c, 'pl.froh'])/sqrt(n),
                mean(froh.stats[froh.stats$covg == c, 'gt.froh']),
                sd(froh.stats[froh.stats$covg == c, 'gt.froh'])/sqrt(n),
                mean(froh.stats[froh.stats$covg == c, 'plink.froh']),
                sd(froh.stats[froh.stats$covg == c, 'plink.froh'])/sqrt(n))
      OUT <- rbind(OUT, save)
    }
    OUT <- as.data.frame(OUT)
    colnames(OUT) <- c('covg','mean.true','se.true','mean.pl','se.pl','mean.gt','se.gt','mean.plink','se.plink')
    # write.csv(OUT, '../manuscript/r_scripts_AMH/tables_output/covg_fROH_means_SEs_simulated.csv')

    k <- 1
    while(k == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_simulated_coverage_vs_fROH_indivlines_n',n,'.pdf'), width = 10, height = 4)
      figure.ct <- figure.ct + 1
      par(mfrow=c(1,3))
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$gt.froh[1]), col = alpha('black', alph))
        lines(c(2:6), c(temp$gt.froh), col = alpha(gt.col, alph))
        points(c(2:6), c(temp$gt.froh), pch = 16, col = alpha(gt.col, pt.alph))
        # points(1, temp$true.froh[1], pch = 16, col = 'black')
        lines(c(2:6), OUT$mean.gt, col = gt.col, lwd = lwd)
        points(c(2:6), OUT$mean.gt, col = gt.col, pch = 23, bg = 'white', cex = 1.5)
      }
      ## PL
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likelihoods - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$pl.froh[1]), col = alpha('black', alph))
        lines(c(2:6), c(temp$pl.froh), col = alpha(pl.col, alph))
        points(c(2:6), c(temp$pl.froh), pch = 16, col = alpha(pl.col, pt.alph))
        # points(1, temp$true.froh[1], pch = 16, col = 'black')
        lines(c(2:6), OUT$mean.pl, col = pl.col, lwd = lwd)
        points(c(2:6), OUT$mean.pl, col = pl.col, pch = 23, bg = 'white', cex = 1.5)
      }
      ## PLINK
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'PLINK', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$plink.froh[1]), col = alpha('black', alph))
        lines(c(2:6), c(temp$plink.froh), col = alpha(plink.col, alph))
        points(c(2:6), c(temp$plink.froh), pch = 16, col = alpha(plink.col, pt.alph))
        # points(1, temp$true.froh[1], pch = 16, col = 'black')
        lines(c(2:6), OUT$mean.plink, col = plink.col, lwd = lwd)
        points(c(2:6), OUT$mean.plink, col = plink.col, pch = 23, bg = 'white', cex = 1.5)
      }
      dev.off()

      ### background individual lines + means and 95% CIs
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_simulated_coverage_vs_fROH_indivlines_n',n,'_means_95CIs.pdf'), width = 10, height = 4)
      figure.ct <- figure.ct + 1
      par(mfrow=c(1,3))

      ## GT
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2, col = true.col)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$gt.froh[1]), col = alpha('black', bg.alph))
        lines(c(2:6), c(temp$gt.froh), col = alpha(gt.col, bg.alph))
      }

      # lines(c(1,2), c(OUT$mean.true[1], OUT$mean.gt[1]), lwd = lwd)
      # points(1, OUT$mean.true[1], pch = 16, col = 'black', cex = pt.size)
      # arrows(x0 = 1, x1 = 1, y0 = c(OUT$mean.true[1] - OUT$se.true[1]*1.96),
      #        y1 = c(OUT$mean.true[1] + OUT$se.true[1]*1.96),
      #        lwd = lwd, col = 'black', code=3, angle=90, length= err.wid)

      lines(c(2:6), c(OUT$mean.gt), col = gt.col, lwd = lwd)
      points(c(2:6), OUT$mean.gt, pch = 16, col = gt.col, cex = pt.size)
      arrows(x0 = c(2:6), x1 = c(2:6), y0 = c(OUT$mean.gt - OUT$se.gt*1.96),
             y1 = c(OUT$mean.gt + OUT$se.gt*1.96),
             lwd = lwd, col = gt.col, code=3, angle=90, length= err.wid)

      ## PL
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likelihoods - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2, col = true.col)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$pl.froh[1]), col = alpha('black', bg.alph))
        lines(c(2:6), c(temp$pl.froh), col = alpha(pl.col, bg.alph))
      }

      # lines(c(1,2), c(OUT$mean.true[1], OUT$mean.pl[1]), lwd = lwd)
      # points(1, OUT$mean.true[1], pch = 16, col = 'black', cex = pt.size)
      # arrows(x0 = 1, x1 = 1, y0 = c(OUT$mean.true[1] - OUT$se.true[1]*1.96),
      #        y1 = c(OUT$mean.true[1] + OUT$se.true[1]*1.96),
      #        lwd = lwd, col = 'black', code=3, angle=90, length= err.wid)

      lines(c(2:6), c(OUT$mean.pl), col = pl.col, lwd = lwd)
      points(c(2:6), OUT$mean.pl, pch = 16, col = pl.col, cex = pt.size)
      arrows(x0 = c(2:6), x1 = c(2:6), y0 = c(OUT$mean.pl - OUT$se.pl*1.96),
             y1 = c(OUT$mean.pl + OUT$se.pl*1.96),
             lwd = lwd, col = pl.col, code=3, angle=90, length= err.wid)

      ## PLINK
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likelihoods - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2, col = true.col)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$plink.froh[1]), col = alpha('black', bg.alph))
        lines(c(2:6), c(temp$plink.froh), col = alpha(plink.col, bg.alph))
      }

      # lines(c(1,2), c(OUT$mean.true[1], OUT$mean.plink[1]), lwd = lwd)
      # points(1, OUT$mean.true[1], pch = 16, col = 'black', cex = pt.size)
      # arrows(x0 = 1, x1 = 1, y0 = c(OUT$mean.true[1] - OUT$se.true[1]*1.96),
      #        y1 = c(OUT$mean.true[1] + OUT$se.true[1]*1.96),
      #        lwd = lwd, col = 'black', code=3, angle=90, length= err.wid)

      lines(c(2:6), c(OUT$mean.plink), col = plink.col, lwd = lwd)
      points(c(2:6), OUT$mean.plink, pch = 16, col = plink.col, cex = pt.size)
      arrows(x0 = c(2:6), x1 = c(2:6), y0 = c(OUT$mean.plink - OUT$se.plink*1.96),
             y1 = c(OUT$mean.plink + OUT$se.plink*1.96),
             lwd = lwd, col = plink.col, code=3, angle=90, length= err.wid)


      dev.off()

      ### background individual lines + means and 83% CIs
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_simulated_coverage_vs_fROH_indivlines_n',n,'_means_83CIs.pdf'), width = 10, height = 4)
      figure.ct <- figure.ct + 1
      par(mfrow=c(1,3))

      ## GT
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2, col = true.col)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$gt.froh[1]), col = alpha('black', bg.alph))
        lines(c(2:6), c(temp$gt.froh), col = alpha(gt.col, bg.alph))
      }

      # lines(c(1,2), c(OUT$mean.true[1], OUT$mean.gt[1]), lwd = lwd)
      # points(1, OUT$mean.true[1], pch = 16, col = 'black', cex = pt.size)
      # arrows(x0 = 1, x1 = 1, y0 = c(OUT$mean.true[1] - OUT$se.true[1]*1.37),
      #        y1 = c(OUT$mean.true[1] + OUT$se.true[1]*1.37),
      #        lwd = lwd, col = 'black', code=3, angle=90, length= err.wid)

      lines(c(2:6), c(OUT$mean.gt), col = gt.col, lwd = lwd)
      points(c(2:6), OUT$mean.gt, pch = 16, col = gt.col, cex = pt.size)
      arrows(x0 = c(2:6), x1 = c(2:6), y0 = c(OUT$mean.gt - OUT$se.gt*1.37),
             y1 = c(OUT$mean.gt + OUT$se.gt*1.37),
             lwd = lwd, col = gt.col, code=3, angle=90, length= err.wid)

      ## PL
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likelihoods - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2, col = true.col)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$pl.froh[1]), col = alpha('black', bg.alph))
        lines(c(2:6), c(temp$pl.froh), col = alpha(pl.col, bg.alph))
      }

      # lines(c(1,2), c(OUT$mean.true[1], OUT$mean.pl[1]), lwd = lwd)
      # points(1, OUT$mean.true[1], pch = 16, col = 'black', cex = pt.size)
      # arrows(x0 = 1, x1 = 1, y0 = c(OUT$mean.true[1] - OUT$se.true[1]*1.37),
      #        y1 = c(OUT$mean.true[1] + OUT$se.true[1]*1.37),
      #        lwd = lwd, col = 'black', code=3, angle=90, length= err.wid)

      lines(c(2:6), c(OUT$mean.pl), col = pl.col, lwd = lwd)
      points(c(2:6), OUT$mean.pl, pch = 16, col = pl.col, cex = pt.size)
      arrows(x0 = c(2:6), x1 = c(2:6), y0 = c(OUT$mean.pl - OUT$se.pl*1.37),
             y1 = c(OUT$mean.pl + OUT$se.pl*1.37),
             lwd = lwd, col = pl.col, code=3, angle=90, length= err.wid)

      ## PLINK
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likelihoods - simulated', cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = txt.size)
      abline(h = OUT$mean.true[1], lty = 2, col = true.col)
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        # lines(c(1:2), c(temp$true.froh[1], temp$plink.froh[1]), col = alpha('black', bg.alph))
        lines(c(2:6), c(temp$plink.froh), col = alpha(plink.col, bg.alph))
      }

      # lines(c(1,2), c(OUT$mean.true[1], OUT$mean.plink[1]), lwd = lwd)
      # points(1, OUT$mean.true[1], pch = 16, col = 'black', cex = pt.size)
      # arrows(x0 = 1, x1 = 1, y0 = c(OUT$mean.true[1] - OUT$se.true[1]*1.37),
      #        y1 = c(OUT$mean.true[1] + OUT$se.true[1]*1.37),
      #        lwd = lwd, col = 'black', code=3, angle=90, length= err.wid)

      lines(c(2:6), c(OUT$mean.plink), col = plink.col, lwd = lwd)
      points(c(2:6), OUT$mean.plink, pch = 16, col = plink.col, cex = pt.size)
      arrows(x0 = c(2:6), x1 = c(2:6), y0 = c(OUT$mean.plink - OUT$se.plink*1.37),
             y1 = c(OUT$mean.plink + OUT$se.plink*1.37),
             lwd = lwd, col = plink.col, code=3, angle=90, length= err.wid)


      dev.off()

      k <- k+1
    }

    ## RELATIVE TO MEAN AT EACH COVERAGE LEVEL
    alph <- 0.5  ## 1 for n = 30, 0.5 for n = 100

    txt.size <- 1.25
    ymin <- -0.5
    ymax <- 0.5
    xmin <- 0.85
    xmax <- 6.15

    k <- 1
    while(k == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_simulated_coverage_vs_fROH_indivlines_n',n,'_REL_TO_MEAN.pdf'), width = 10, height = 4)
      figure.ct <- figure.ct + 1
      par(mfrow=c(1,3))

      ## GT
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'gt.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$gt.froh <- temp$gt.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$gt.froh[1]), col = alpha('black', alph))
        lines(c(2:6), c(temp$gt.froh), col = alpha(gt.col, alph))
        points(c(2:6), c(temp$gt.froh), pch = 16, col = gt.col)

        ## only plotting 15X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$gt.froh[3]), col = alpha('black', alph))
        # lines(c(2:4), c(temp$gt.froh[3:5]), col = alpha(gt.col, alph))
        # points(c(2:4), c(temp$gt.froh[3:5]), pch = 16, col = gt.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      ## PL
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likeilhoods - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'pl.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$pl.froh <- temp$pl.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$pl.froh[1]), col = alpha('black', alph))
        lines(c(2:6), c(temp$pl.froh), col = alpha(pl.col, alph))
        points(c(2:6), c(temp$pl.froh), pch = 16, col = pl.col)

        ## only plotting 15X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$pl.froh[3]), col = alpha('black', alph))
        # lines(c(2:4), c(temp$pl.froh[3:5]), col = alpha(pl.col, alph))
        # points(c(2:4), c(temp$pl.froh[3:5]), pch = 16, col = pl.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      ## PLINK
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'PLINK - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'plink.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$plink.froh <- temp$plink.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$plink.froh[1]), col = alpha('black', alph))
        lines(c(2:6), c(temp$plink.froh), col = alpha(plink.col, alph))
        points(c(2:6), c(temp$plink.froh), pch = 16, col = plink.col)

        ## only plotting 15X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$plink.froh[3]), col = alpha('black', alph))
        # lines(c(2:4), c(temp$plink.froh[3:5]), col = alpha(plink.col, alph))
        # points(c(2:4), c(temp$plink.froh[3:5]), pch = 16, col = plink.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }


      ## GT
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'gt.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$gt.froh <- temp$gt.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$gt.froh[1]), col = alpha('black', alph))
        # lines(c(2:6), c(temp$gt.froh), col = alpha(gt.col, alph))
        # points(c(2:6), c(temp$gt.froh), pch = 16, col = gt.col)

        ## only plotting 15X - 50X
        lines(c(1,4), c((temp$true.froh[1] - MEANS[1]), temp$gt.froh[3]), col = alpha('black', alph))
        lines(c(4:6), c(temp$gt.froh[3:5]), col = alpha(gt.col, alph))
        points(c(4:6), c(temp$gt.froh[3:5]), pch = 16, col = gt.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      ## PL
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likeilhoods - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'pl.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$pl.froh <- temp$pl.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$pl.froh[1]), col = alpha('black', alph))
        # lines(c(2:6), c(temp$pl.froh), col = alpha(pl.col, alph))
        # points(c(2:6), c(temp$pl.froh), pch = 16, col = pl.col)

        ## only plotting 15X - 50X
        lines(c(1,4), c((temp$true.froh[1] - MEANS[1]), temp$pl.froh[3]), col = alpha('black', alph))
        lines(c(4:6), c(temp$pl.froh[3:5]), col = alpha(pl.col, alph))
        points(c(4:6), c(temp$pl.froh[3:5]), pch = 16, col = pl.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      ## PLINK
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'PLINK - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'plink.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$plink.froh <- temp$plink.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$plink.froh[1]), col = alpha('black', alph))
        # lines(c(2:6), c(temp$plink.froh), col = alpha(plink.col, alph))
        # points(c(2:6), c(temp$plink.froh), pch = 16, col = plink.col)

        ## only plotting 15X - 50X
        lines(c(1,4), c((temp$true.froh[1] - MEANS[1]), temp$plink.froh[3]), col = alpha('black', alph))
        lines(c(4:6), c(temp$plink.froh[3:5]), col = alpha(plink.col, alph))
        points(c(4:6), c(temp$plink.froh[3:5]), pch = 16, col = plink.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      ## GT
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Genotypes - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'gt.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$gt.froh <- temp$gt.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$gt.froh[1]), col = alpha('black', alph))
        # lines(c(2:6), c(temp$gt.froh), col = alpha(gt.col, alph))
        # points(c(2:6), c(temp$gt.froh), pch = 16, col = gt.col)

        ## only plotting 15X - 50X
        lines(c(1,6), c((temp$true.froh[1] - MEANS[1]), temp$gt.froh[5]), col = alpha('black', alph))
        points(c(6), c(temp$gt.froh[5]), pch = 16, col = gt.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      ## PL
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'Likeilhoods - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'pl.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$pl.froh <- temp$pl.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$pl.froh[1]), col = alpha('black', alph))
        # lines(c(2:6), c(temp$pl.froh), col = alpha(pl.col, alph))
        # points(c(2:6), c(temp$pl.froh), pch = 16, col = pl.col)

        ## only plotting 15X - 50X
        lines(c(1,6), c((temp$true.froh[1] - MEANS[1]), temp$pl.froh[5]), col = alpha('black', alph))
        points(c(6), c(temp$pl.froh[5]), pch = 16, col = pl.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      ## PLINK
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ylab = 'f(ROH)',
           xaxt = 'n', xlab = 'Coverage', main = 'PLINK - simulated', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = txt.size)
      # axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c('0.0','0.1','0.2','0.3','0.4','0.5'), cex.axis = txt.size)
      MEANS <- mean(froh.stats[froh.stats$covg == 5, 'true.froh'])
      for(c in sort(unique(froh.stats$covg))){
        MEANS <- c(MEANS, mean(froh.stats[froh.stats$covg == c, 'plink.froh']))
      }
      abline(h = 0, lty = 2, col = 'lightgray')
      for(i in unique(froh.stats$id)){
        temp <- froh.stats[froh.stats$id == i,]
        temp <- temp[order(temp$covg),]
        temp$plink.froh <- temp$plink.froh - MEANS[c(2:6)]

        ## plotting 5X - 50X
        # lines(c(1:2), c((temp$true.froh[1] - MEANS[1]), temp$plink.froh[1]), col = alpha('black', alph))
        # lines(c(2:6), c(temp$plink.froh), col = alpha(plink.col, alph))
        # points(c(2:6), c(temp$plink.froh), pch = 16, col = plink.col)

        ## only plotting 15X - 50X
        lines(c(1,6), c((temp$true.froh[1] - MEANS[1]), temp$plink.froh[5]), col = alpha('black', alph))
        points(c(6), c(temp$plink.froh[5]), pch = 16, col = plink.col)

        points(1, (temp$true.froh[1] - MEANS[1]), pch = 16, col = 'black')
      }

      dev.off()
      k <- k+1
    }


    ##### 93. Looking at BCFtools ROH quality measure #####
    gt.overlap$prop.true <- gt.overlap$true.len/gt.overlap$len
    pl.overlap$prop.true <- pl.overlap$true.len/pl.overlap$len

    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_roh_quality_measure_test.pdf'), width = 12, height = 5)
    figure.ct <- figure.ct + 1
    par(mfrow = c(1,2))
    for(c in sort(unique(gt.overlap$covg))){
      sub <- gt.overlap[gt.overlap$covg == c,]
      plot(sub$prop.true, sub$roh.call.quality, pch = 16, cex = 0.5, col = alpha('blue', 0.2), main = paste0('GT - ',c),
           xlab = 'Proportion of called ROH that overlaps true ROH(s)', ylab = 'ROH call quality score')
      plot(sub$len, sub$prop.true, pch = 16, col = alpha('blue', 0.2), cex = 0.5, main = paste0('GT - ',c),
           xlab = 'Called ROH length', ylab = 'Proportion of called ROH that overlaps true ROH(s)')
      pts <- seq(1,4e6,1e3)
      pts <- cbind(pts, (1e5/pts))
      lines(pts[,1], pts[,2], lty = 2)
    }

    for(c in sort(unique(pl.overlap$covg))){
      sub <- pl.overlap[pl.overlap$covg == c,]
      plot(sub$prop.true, sub$roh.call.quality, pch = 16, cex = 0.5, col = alpha('blue', 0.2), main = paste0('PL - ',c),
           xlab = 'Proportion of called ROH that overlaps true ROH(s)', ylab = 'ROH call quality score')
      plot(sub$len, sub$prop.true, pch = 16, col = alpha('blue', 0.2), cex = 0.5, main = paste0('PL - ',c),
           xlab = 'Called ROH length', ylab = 'Proportion of called ROH that overlaps true ROH(s)')
      pts <- seq(1,4e6,1e3)
      pts <- cbind(pts, (1e5/pts))
      lines(pts[,1], pts[,2], lty = 2)
    }
    dev.off()

    ##### 94A. Plotting f(ROH) by length bins, individual lines #####
    b1 <- 500e3
    b2 <- 1e6
    b3 <- 2e6

    ymin <- 0
    ymax <- 1
    alph <- 0.3

    k <- 1
    while(k == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_fROH_by_length_bins_indivlines_',n,'.pdf'), width = 16, height = 12)
      figure.ct <- figure.ct + 1
      par(mfrow = c(3,4))

      ## GT
      OUT <- NULL
      for(i in unique(bcf.gt.res$id)){
        for(c in unique(bcf.gt.res$covg)){
          temp <- bcf.gt.res[bcf.gt.res$id == i & bcf.gt.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      gt.frohs <- as.data.frame(OUT)
      colnames(gt.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(gt.frohs)){
        gt.frohs[,c] <- as.numeric(gt.frohs[,c])
      }
      write.csv(gt.frohs, paste0('3_methods_results/',d,'_GT_frohs_by_bins.csv'), row.names = FALSE)

      ## GT plot
      # par(mfrow = c(2,4))
      ## bin1
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotypes only: Short ROHs', xlab = 'Coverage', ylab = 'f(ROH)',
           cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len, temp$bin1.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin1.froh, pch = 16, col = gt.col)
        lines(c(2:6), temp$bin1.froh, col = alpha(gt.col, alph))
      }
      ## bin 2
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotypes only: Intermediate ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len, temp$bin2.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin2.froh, pch = 16, col = gt.col)
        lines(c(2:6), temp$bin2.froh, col = alpha(gt.col, alph))
      }
      ## bin 3
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotypes only: Long ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len, temp$bin3.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin3.froh, pch = 16, col = gt.col)
        lines(c(2:6), temp$bin3.froh, col = alpha(gt.col, alph))
      }
      ## bin 4
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotypes only: Longer ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b3, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length > b3, 'length'])/chrom.len, temp$bin4.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin4.froh, pch = 16, col = gt.col)
        lines(c(2:6), temp$bin4.froh, col = alpha(gt.col, alph))
      }

      ## PL
      OUT <- NULL
      for(i in unique(bcf.pl.res$id)){
        for(c in unique(bcf.pl.res$covg)){
          temp <- bcf.pl.res[bcf.pl.res$id == i & bcf.pl.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      pl.frohs <- as.data.frame(OUT)
      colnames(pl.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(pl.frohs)){
        pl.frohs[,c] <- as.numeric(pl.frohs[,c])
      }
      write.csv(pl.frohs, paste0('3_methods_results/',d,'_PL_frohs_by_bins.csv'), row.names = FALSE)

      ## PL plot
      ## bin1
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotype likelihoods: Short ROHs', xlab = 'Coverage', ylab = 'f(ROH)',
           cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len, temp$bin1.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin1.froh, pch = 16, col = pl.col)
        lines(c(2:6), temp$bin1.froh, col = alpha(pl.col, alph))
      }
      ## bin 2
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotype likelihoods: Intermediate ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len, temp$bin2.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin2.froh, pch = 16, col = pl.col)
        lines(c(2:6), temp$bin2.froh, col = alpha(pl.col, alph))
      }
      ## bin 3
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotype likelihoods: Long ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len, temp$bin3.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin3.froh, pch = 16, col = pl.col)
        lines(c(2:6), temp$bin3.froh, col = alpha(pl.col, alph))
      }
      ## bin 4
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotype likelihoods: Longer ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b3, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length > b3, 'length'])/chrom.len, temp$bin4.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin4.froh, pch = 16, col = pl.col)
        lines(c(2:6), temp$bin4.froh, col = alpha(pl.col, alph))
      }

      ## PLINK
      OUT <- NULL
      for(i in unique(plink.res$id)){
        for(c in unique(plink.res$covg)){
          temp <- plink.res[plink.res$id == i & plink.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      plink.frohs <- as.data.frame(OUT)
      colnames(plink.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(plink.frohs)){
        plink.frohs[,c] <- as.numeric(plink.frohs[,c])
      }
      write.csv(plink.frohs, paste0('3_methods_results/',d,'_PLINK_frohs_by_bins.csv'), row.names = FALSE)

      ## PLINK plot
      ## bin1
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'PLINK: Short ROHs', xlab = 'Coverage', ylab = 'f(ROH)',
           cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len, temp$bin1.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin1.froh, pch = 16, col = plink.col)
        lines(c(2:6), temp$bin1.froh, col = alpha(plink.col, alph))
      }
      ## bin 2
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'PLINK: Intermediate ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len, temp$bin2.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin2.froh, pch = 16, col = plink.col)
        lines(c(2:6), temp$bin2.froh, col = alpha(plink.col, alph))
      }
      ## bin 3
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'PLINK: Long ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len, temp$bin3.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin3.froh, pch = 16, col = plink.col)
        lines(c(2:6), temp$bin3.froh, col = alpha(plink.col, alph))
      }
      ## bin 4
      plot(0,0, xlim = c(1,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'PLINK: Longer ROHs', xlab = 'Coverage', ylab = 'f(ROH)', cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(1,2,3,4,5,6), labels = c('True','5X','10X','15X','30X','50X'), cex.axis = 1.25)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        points(1, sum(true.rohs[true.rohs$id == i & true.rohs$length >= b3, 'length'])/chrom.len, pch = 16, col = 'black')
        lines(x = c(1,2), y = c(sum(true.rohs[true.rohs$id == i & true.rohs$length > b3, 'length'])/chrom.len, temp$bin4.froh[1]), col = alpha('black', alph))
        points(c(2:6), temp$bin4.froh, pch = 16, col = plink.col)
        lines(c(2:6), temp$bin4.froh, col = alpha(plink.col, alph))
      }

      k <- k+1
      dev.off()
    }


    ##### >>> 94B. Plotting f(ROH) by length bins, individual lines - relative values (formatting for SI) #####
    b1 <- 500e3
    b2 <- 1e6
    b3 <- 2e6
    
    if(d == 'decline'){
      ymin <- -0.55
      ymax <- 0.45
    }
    if(d == 'bottle'){
      ymin <- -0.4
      ymax <- 0.45
    }
    if(d == 'small'){
      ymin <- -0.6
      ymax <- 0.35
    }
    if(d == 'large-1000'){
      ymin <- -0.25
      ymax <- 0.15
    }

    ln.alph <- 0.3
    pt.alph <- 0.5
    txt.size <- 2

    k <- 1
    while(k == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_fROH_by_length_bins_indivlines_',n,'_relative_vals.pdf'), width = 12, height = 10)
      figure.ct <- figure.ct + 1
      par(mfrow = c(3,4), mar = c(2.5, 2.5, 0.6, 0.6))

      ## GT
      OUT <- NULL
      for(i in unique(bcf.gt.res$id)){
        for(c in unique(bcf.gt.res$covg)){
          temp <- bcf.gt.res[bcf.gt.res$id == i & bcf.gt.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
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
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '',
           cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len
        points(c(2:6), temp$bin1.froh - true.val, pch = 16, col = alpha(gt.col, pt.alph))
        lines(c(2:6), temp$bin1.froh - true.val, col = alpha(gt.col, ln.alph))
      }
      ## bin 2
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len
        points(c(2:6), temp$bin2.froh - true.val, pch = 16, col = alpha(gt.col, pt.alph))
        lines(c(2:6), temp$bin2.froh - true.val, col = alpha(gt.col, ln.alph))
      }
      ## bin 3
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len
        points(c(2:6), temp$bin3.froh - true.val, pch = 16, col = alpha(gt.col, pt.alph))
        lines(c(2:6), temp$bin3.froh - true.val, col = alpha(gt.col, ln.alph))
      }
      ## bin 4
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(gt.frohs$id)){
        temp <- gt.frohs[gt.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length > b3, 'length'])/chrom.len
        points(c(2:6), temp$bin4.froh - true.val, pch = 16, col = alpha(gt.col, pt.alph))
        lines(c(2:6), temp$bin4.froh - true.val, col = alpha(gt.col, ln.alph))
      }

      ## PL
      OUT <- NULL
      for(i in unique(bcf.pl.res$id)){
        for(c in unique(bcf.pl.res$covg)){
          temp <- bcf.pl.res[bcf.pl.res$id == i & bcf.pl.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
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
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '',
           cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len
        points(c(2:6), temp$bin1.froh - true.val, pch = 16, col = alpha(pl.col, pt.alph))
        lines(c(2:6), temp$bin1.froh - true.val, col = alpha(pl.col, ln.alph))
      }
      ## bin 2
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len
        points(c(2:6), temp$bin2.froh - true.val, pch = 16, col = alpha(pl.col, pt.alph))
        lines(c(2:6), temp$bin2.froh - true.val, col = alpha(pl.col, ln.alph))
      }
      ## bin 3
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len
        points(c(2:6), temp$bin3.froh - true.val, pch = 16, col = alpha(pl.col, pt.alph))
        lines(c(2:6), temp$bin3.froh - true.val, col = alpha(pl.col, ln.alph))
      }
      ## bin 4
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      # axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(pl.frohs$id)){
        temp <- pl.frohs[pl.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length > b3, 'length'])/chrom.len
        points(c(2:6), temp$bin4.froh - true.val, pch = 16, col = alpha(pl.col, pt.alph))
        lines(c(2:6), temp$bin4.froh - true.val, col = alpha(pl.col, ln.alph))
      }

      ## PLINK
      OUT <- NULL
      for(i in unique(plink.res$id)){
        for(c in unique(plink.res$covg)){
          temp <- plink.res[plink.res$id == i & plink.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
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
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '',
           cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size, line = 0.5, lwd = 0)
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length < b1, 'length'])/chrom.len
        points(c(2:6), temp$bin1.froh - true.val, pch = 16, col = alpha(plink.col, pt.alph))
        lines(c(2:6), temp$bin1.froh - true.val, col = alpha(plink.col, ln.alph))
      }
      ## bin 2
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size, line = 0.5, lwd = 0)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length >= b1 & true.rohs$length < b2, 'length'])/chrom.len
        points(c(2:6), temp$bin2.froh - true.val, pch = 16, col = alpha(plink.col, pt.alph))
        lines(c(2:6), temp$bin2.froh - true.val, col = alpha(plink.col, ln.alph))
      }
      ## bin 3
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size, line = 0.5, lwd = 0)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length >= b2 & true.rohs$length < b3, 'length'])/chrom.len
        points(c(2:6), temp$bin3.froh - true.val, pch = 16, col = alpha(plink.col, pt.alph))
        lines(c(2:6), temp$bin3.froh - true.val, col = alpha(plink.col, ln.alph))
      }
      ## bin 4
      plot(0,0, xlim = c(2,6), ylim = c(ymin, ymax),
           xaxt = 'n', main = '', xlab = '', ylab = '', cex.axis = txt.size, cex.lab = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('','','','',''), cex.axis = txt.size)
      axis(1, at = c(2,3,4,5,6), labels = c('5X','10X','15X','30X','50X'), cex.axis = txt.size, line = 0.5, lwd = 0)
      # axis(2, at = c(-0.2, -0.1, 0.0, 0.1, 0.2), labels = c('','','','',''))
      abline(h = 0, lty = 2)
      x <- 1
      for(i in unique(plink.frohs$id)){
        temp <- plink.frohs[plink.frohs$id == i,]
        temp <- temp[order(temp$covg),]
        true.val <- sum(true.rohs[true.rohs$id == i & true.rohs$length > b3, 'length'])/chrom.len
        points(c(2:6), temp$bin4.froh - true.val, pch = 16, col = alpha(plink.col, pt.alph))
        lines(c(2:6), temp$bin4.froh - true.val, col = alpha(plink.col, ln.alph))
      }

      k <- k+1
      dev.off()
    }

    ### for each individual, calculate how often the 5X estimate is > the 10X estimate
    OUT <- NULL
    for(i in unique(pl.frohs$id)){
      sub.pl <- pl.frohs[pl.frohs$id == i,]
      sub.gt <- gt.frohs[gt.frohs$id == i,]
      sub.plink <- plink.frohs[plink.frohs$id == i,]
      save <- c(i, sub.pl[sub.pl$covg == 5, 'bin1.froh'] - sub.pl[sub.pl$covg == 10, 'bin1.froh'],
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
    OUT <- OUT[,-1]
    cats <- c('PL short','PL intermediate','PL long','PL very long',
              'GT short','GT intermediate','GT long','GT very long',
              'PLINK short','PLINK intermediate','PLINK long','PLINK very long')
    print('Proportion of 5X f(ROH) estimates > 10X f(ROH) estimates')
    for(c in 1:ncol(OUT)){
      print(paste0(cats[c],': ',nrow(OUT[OUT[,c] > 0,])/nrow(OUT)))
    }
    print('Proportion of 5X f(ROH) estimates >= 10X f(ROH) estimates')
    for(c in 1:ncol(OUT)){
      print(paste0(cats[c],': ',nrow(OUT[OUT[,c] >= 0,])/nrow(OUT)))
    }

    ##### >>> 94C. Plotting f(ROH) by length bins, 1 plot per method, lines = different coverages #####
    b1 <- 500e3
    b2 <- 1e6
    b3 <- 2e6

    n <- length(unique(froh.stats$id))

    k <- 1
    while(k == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_fROH_by_length_bins_indivlines_',n,'_relative_vals_all3methods.pdf'), width = 6, height = 5)
      figure.ct <- figure.ct + 1

      ## GT
      OUT <- NULL
      for(i in unique(bcf.gt.res$id)){
        for(c in unique(bcf.gt.res$covg)){
          temp <- bcf.gt.res[bcf.gt.res$id == i & bcf.gt.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
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
      for(i in unique(bcf.pl.res$id)){
        for(c in unique(bcf.pl.res$covg)){
          temp <- bcf.pl.res[bcf.pl.res$id == i & bcf.pl.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
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
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      plink.frohs <- as.data.frame(OUT)
      colnames(plink.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(plink.frohs)){
        plink.frohs[,c] <- as.numeric(plink.frohs[,c])
      }

      ## calculate true f(ROH) for bins
      OUT <- NULL
      for(i in unique(true.rohs$id)){
        temp <- true.rohs[true.rohs$id == i,]
        true.total <- sum(temp$length)/chrom.len
        true.b1 <- sum(temp[temp$length < b1, 'length'])/chrom.len
        true.b2 <- sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len
        true.b3 <- sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len
        true.b4 <- sum(temp[temp$length >= b3, 'length'])/chrom.len
        save <- c(i, true.total, true.b1, true.b2, true.b3, true.b4)
        OUT <- rbind(OUT, save)
      }
      true.roh.bins <- as.data.frame(OUT)
      colnames(true.roh.bins) <- c('id','true.total','true.b1','true.b2','true.b3','true.b4')
      for(c in 1:ncol(true.roh.bins)){
        true.roh.bins[,c] <- as.numeric(true.roh.bins[,c])
      }

      ## Combined plot
      ymin <- -0.5
      ymax <- 0.5
      ln.alph <- 0.5
      pt.alph <- 1
      diff <- 0.15
      xmin <- 1.75
      xmax <- 5.25
      orig.xs <- c(2:5)

      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotypes only', xlab = 'Coverage', 
           ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])),
           cex.axis = 1.25, cex.lab = 1.25)
      axis(1, at = c(2,3,4,5), labels = c('Short','Intermediate','Long','Very long'), cex.axis = 1.25)
      abline(h = 0, lty = 2)

      ## GT
      col <- 1
      for(c in sort(unique(gt.frohs$covg))){
        temp <- gt.frohs[gt.frohs$covg == c,]
        temp <- merge(temp, true.roh.bins, by = 'id')
        temp$b1.diff <- temp$bin1.froh - temp$true.b1
        temp$b2.diff <- temp$bin2.froh - temp$true.b2
        temp$b3.diff <- temp$bin3.froh - temp$true.b3
        temp$b4.diff <- temp$bin4.froh - temp$true.b4

        xs <- orig.xs
        columns <- c(18, 19, 20, 21)
        lines(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(gt.cols[col], ln.alph))
        points(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(gt.cols[col], pt.alph), pch = 16)
        for(l in c(1:4)){
          column <- columns[l]
          arrows(x0 = xs[l], x1 = xs[l], y0 = (mean(temp[,column], na.rm = TRUE) - sd(temp[,column], na.rm = TRUE)),
                 y1 = (mean(temp[,column], na.rm = TRUE) + sd(temp[,column], na.rm = TRUE)),
                 lwd = 2, col = alpha(gt.cols[col], pt.alph), code=3, angle=90, length=0.1)
        }
        col <- col+1
      }

      ## PL
      col <- 1
      for(c in sort(unique(pl.frohs$covg))){
        temp <- pl.frohs[pl.frohs$covg == c,]
        temp <- merge(temp, true.roh.bins, by = 'id')
        temp$b1.diff <- temp$bin1.froh - temp$true.b1
        temp$b2.diff <- temp$bin2.froh - temp$true.b2
        temp$b3.diff <- temp$bin3.froh - temp$true.b3
        temp$b4.diff <- temp$bin4.froh - temp$true.b4

        xs <- orig.xs - diff
        lines(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(pl.cols[col], ln.alph))
        points(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(pl.cols[col], pt.alph), pch = 16)
        for(l in c(1:4)){
          column <- columns[l]
          arrows(x0 = xs[l], x1 = xs[l], y0 = (mean(temp[,column], na.rm = TRUE) - sd(temp[,column], na.rm = TRUE)),
                 y1 = (mean(temp[,column], na.rm = TRUE) + sd(temp[,column], na.rm = TRUE)),
                 lwd = 2, col = alpha(pl.cols[col], pt.alph), code=3, angle=90, length=0.1)
        }
        col <- col+1

        col <- col+1
      }

      ## PLINK
      col <- 1
      for(c in sort(unique(plink.frohs$covg))){
        temp <- plink.frohs[plink.frohs$covg == c,]
        temp <- merge(temp, true.roh.bins, by = 'id')
        temp$b1.diff <- temp$bin1.froh - temp$true.b1
        temp$b2.diff <- temp$bin2.froh - temp$true.b2
        temp$b3.diff <- temp$bin3.froh - temp$true.b3
        temp$b4.diff <- temp$bin4.froh - temp$true.b4

        xs <- orig.xs + diff
        lines(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(plink.cols[col], ln.alph))
        points(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(plink.cols[col], pt.alph), pch = 16)
        for(l in c(1:4)){
          column <- columns[l]
          arrows(x0 = xs[l], x1 = xs[l], y0 = (mean(temp[,column], na.rm = TRUE) - sd(temp[,column], na.rm = TRUE)),
                 y1 = (mean(temp[,column], na.rm = TRUE) + sd(temp[,column], na.rm = TRUE)),
                 lwd = 2, col = alpha(plink.cols[col], pt.alph), code=3, angle=90, length=0.1)
        }
        col <- col+1
      }

      dev.off()
      k <- k+1
    }


    ### Separate plot per method
    ## --- 95% CIs are comically small, may need to calculate a different way?
    k <- 1
    while(k == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_fROH_by_length_bins_indivlines_',n,'_relative_vals_separate_method_plots_95CIs.pdf'), 
          width = 8, height = 5.5)
      par(mar = c(5.1, 5.1, 4.1, 2.1))
      figure.ct <- figure.ct + 1

      ymin <- -0.55
      ymax <- 0.55
      ln.alph <- 0.5
      pt.alph <- 1
      diff <- 0.15
      xmin <- 1.75
      xmax <- 5.25
      offsets <- c(-0.2, -0.1, 0, 0.1, 0.2)
      orig.xs <- c(2:5)
      text.size <- 1.75
      pt.cex <- 1.25
      lwd <- 2

      ## GT
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Genotypes', xlab = 'ROH length bin', 
           ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])),
           cex.axis = text.size, cex.lab = text.size, yaxt = 'n')
      axis(2, at = c(-0.5, -0.25, 0, 0.25, 0.5), cex.axis = text.size)
      axis(1, at = c(2,3,4,5), labels = c('Short','Intermediate','Long','Very long'), cex.axis = text.size)
      abline(h = 0, lty = 2)

      col <- 1
      for(c in sort(unique(gt.frohs$covg))){
        temp <- gt.frohs[gt.frohs$covg == c,]
        temp <- merge(temp, true.roh.bins, by = 'id')
        temp$b1.diff <- temp$bin1.froh - temp$true.b1
        temp$b2.diff <- temp$bin2.froh - temp$true.b2
        temp$b3.diff <- temp$bin3.froh - temp$true.b3
        temp$b4.diff <- temp$bin4.froh - temp$true.b4

        xs <- orig.xs + offsets[col]
        columns <- c(18, 19, 20, 21)
        lines(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(gt.cols[col], ln.alph), lwd = lwd)
        points(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(gt.cols[col], pt.alph), pch = 16, cex = pt.cex)
        for(l in c(1:4)){
          column <- columns[l]
          ## 95% CIs (inappropriate for large sample sizes)
          # arrows(x0 = xs[l], x1 = xs[l], y0 = (mean(temp[,column], na.rm = TRUE) - (sd(temp[,column], na.rm = TRUE)/10*1.96)),
          #        y1 = (mean(temp[,column], na.rm = TRUE) + (sd(temp[,column], na.rm = TRUE)/10*1.96)),
          #        lwd = 2, col = alpha(gt.cols[col], pt.alph), code=3, angle=90, length=0.1)

          arrows(x0 = xs[l], x1 = xs[l], y0 = quantile(temp[,column], probs = c(0.025,0.975))[1],
                 y1 = quantile(temp[,column], probs = c(0.025,0.975))[2],
                 lwd = 2, col = alpha(gt.cols[col], pt.alph), code=3, angle=90, length=0.1)
        }
        col <- col+1
      }
      legend('top', legend = c('5X','10X','15X','30X','50X'), col = gt.cols, pch = 16, bty = 'n', cex = text.size, pt.cex = pt.cex, horiz = TRUE, x.intersp = 0.7)

      ## PL
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'Likelihoods', xlab = 'ROH length bin', 
           ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])),
           cex.axis = text.size, cex.lab = text.size, yaxt  = 'n')
      axis(2, at = c(-0.5, -0.25, 0, 0.25, 0.5), cex.axis = text.size)
      axis(1, at = c(2,3,4,5), labels = c('Short','Intermediate','Long','Very long'), cex.axis = text.size)
      abline(h = 0, lty = 2)

      col <- 1
      for(c in sort(unique(pl.frohs$covg))){
        temp <- pl.frohs[pl.frohs$covg == c,]
        temp <- merge(temp, true.roh.bins, by = 'id')
        temp$b1.diff <- temp$bin1.froh - temp$true.b1
        temp$b2.diff <- temp$bin2.froh - temp$true.b2
        temp$b3.diff <- temp$bin3.froh - temp$true.b3
        temp$b4.diff <- temp$bin4.froh - temp$true.b4

        xs <- orig.xs + offsets[col]
        columns <- c(18, 19, 20, 21)
        lines(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(pl.cols[col], ln.alph), lwd = lwd)
        points(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(pl.cols[col], pt.alph), pch = 16, cex = pt.cex)
        for(l in c(1:4)){
          column <- columns[l]
          arrows(x0 = xs[l], x1 = xs[l], y0 = quantile(temp[,column], probs = c(0.025,0.975))[1],
                 y1 = quantile(temp[,column], probs = c(0.025,0.975))[2],
                 lwd = 2, col = alpha(pl.cols[col], pt.alph), code=3, angle=90, length=0.1)
        }
        col <- col+1
      }
      legend('top', legend = c('5X','10X','15X','30X','50X'), col = pl.cols, pch = 16, bty = 'n', cex = text.size, pt.cex = pt.cex, horiz = TRUE, x.intersp = 0.7)

      ## PLINK
      plot(0,0, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
           xaxt = 'n', main = 'PLINK', xlab = 'ROH length bin', 
           ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])),
           cex.axis = text.size, cex.lab = text.size, yaxt = 'n')
      axis(2, at = c(-0.5, -0.25, 0, 0.25, 0.5), cex.axis = text.size)
      axis(1, at = c(2,3,4,5), labels = c('Short','Intermediate','Long','Very long'), cex.axis = text.size)
      abline(h = 0, lty = 2)

      col <- 1
      for(c in sort(unique(plink.frohs$covg))){
        temp <- plink.frohs[plink.frohs$covg == c,]
        temp <- merge(temp, true.roh.bins, by = 'id')
        temp$b1.diff <- temp$bin1.froh - temp$true.b1
        temp$b2.diff <- temp$bin2.froh - temp$true.b2
        temp$b3.diff <- temp$bin3.froh - temp$true.b3
        temp$b4.diff <- temp$bin4.froh - temp$true.b4

        xs <- orig.xs + offsets[col]
        columns <- c(18, 19, 20, 21)
        lines(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(plink.cols[col], ln.alph), lwd = lwd)
        points(xs, c(mean(temp$b1.diff), mean(temp$b2.diff), mean(temp$b3.diff), mean(temp$b4.diff)), col = alpha(plink.cols[col], pt.alph), pch = 16, cex = pt.cex)
        for(l in c(1:4)){
          column <- columns[l]
          arrows(x0 = xs[l], x1 = xs[l], y0 = quantile(temp[,column], probs = c(0.025,0.975))[1],
                 y1 = quantile(temp[,column], probs = c(0.025,0.975))[2],
                 lwd = 2, col = alpha(plink.cols[col], pt.alph), code=3, angle=90, length=0.1)
        }
        col <- col+1
      }
      legend('top', legend = c('5X','10X','15X','30X','50X'), col = plink.cols, pch = 16, bty = 'n', cex = text.size, pt.cex = pt.cex, horiz = TRUE, x.intersp = 0.7)

      dev.off()
      k <- k+1
    }

    ### Create true f(ROH) bin distributions
    # pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_dists.pdf'), width = 4, height = 4)
    # plot(density(true.roh.bins$true.b1), xlab = 'True f(ROH) - bin 1', xlim = c(-0.02, 0.25), xaxt = 'n')
    #   polygon(density(true.roh.bins$true.b1), col = 'grey')
    #   axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    #   axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    #   abline(v = mean(true.roh.bins$true.b1), lty = 2, lwd = 3)
    # plot(density(true.roh.bins$true.b2), xlab = 'True f(ROH) - bin 2', xlim = c(-0.02, 0.25), xaxt = 'n')
    #   polygon(density(true.roh.bins$true.b2), col = 'grey')
    #   axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    #   axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    #   abline(v = mean(true.roh.bins$true.b2), lty = 2, lwd = 3)
    # plot(density(true.roh.bins$true.b3), xlab = 'True f(ROH) - bin 3', xlim = c(-0.02, 0.25), xaxt = 'n')
    #   polygon(density(true.roh.bins$true.b3), col = 'grey')
    #   axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    #   axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    #   abline(v = mean(true.roh.bins$true.b3), lty = 2, lwd = 3)
    # plot(density(true.roh.bins$true.b4), xlab = 'True f(ROH) - bin 4', xlim = c(-0.02, 0.25), xaxt = 'n')
    #   polygon(density(true.roh.bins$true.b4), col = 'grey')
    #   axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    #   axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    #   abline(v = mean(true.roh.bins$true.b4), lty = 2, lwd = 3)
    # dev.off()

    min.bin.col <- ghibli_palette('LaputaLight')[5]
    max.bin.col <- ghibli_palette('LaputaDark')[2]

    bin.pal <- colorRampPalette(c(min.bin.col, max.bin.col))
    bin.cols <- bin.pal(4)
    alph <- 0.6
    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_dists_single_plot.pdf'), width = 8, height = 5)
    figure.ct <- figure.ct + 1
    plot(density(true.roh.bins$true.b3), xlab = substitute(paste('True ',italic('F')[ROH])), col = 'transparent',
         cex.axis = text.size, cex.lab = text.size, zero.line = FALSE)
      polygon(density(true.roh.bins$true.b4), col = alpha(bin.cols[4], alph), border = NA)
      lines(density(true.roh.bins$true.b4), col = bin.cols[4])
      polygon(density(true.roh.bins$true.b3), col = alpha(bin.cols[3], alph), border = NA)
      lines(density(true.roh.bins$true.b3), col = bin.cols[3])
      polygon(density(true.roh.bins$true.b2), col = alpha(bin.cols[2], alph), border = NA)
      lines(density(true.roh.bins$true.b2), col = bin.cols[2])
      polygon(density(true.roh.bins$true.b1), col = alpha(bin.cols[1], alph), border = NA)
      lines(density(true.roh.bins$true.b1), col = bin.cols[1])
      abline(v = mean(true.roh.bins$true.b1), col = bin.cols[1], lty = 2, lwd = 2)
      abline(v = mean(true.roh.bins$true.b2), col = bin.cols[2], lty = 2, lwd = 2)
      abline(v = mean(true.roh.bins$true.b3), col = bin.cols[3], lty = 2, lwd = 2)
      abline(v = mean(true.roh.bins$true.b4), col = bin.cols[4], lty = 2, lwd = 2)
      legend('topright', legend = c('Short','Intermediate','Long','Very long'), fill = c(alpha(bin.cols[1], alph),
                                                                                         alpha(bin.cols[2], alph),
                                                                                         alpha(bin.cols[3], alph),
                                                                                         alpha(bin.cols[4], alph)),
             border = c(bin.cols[1:4]), inset = 0.05, bty = 'n', cex = text.size)
    dev.off()


    # polygon(density(true.roh.bins$true.b1), col = 'grey')
    # axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    # axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    # abline(v = mean(true.roh.bins$true.b1), lty = 2, lwd = 3)
    # plot(density(true.roh.bins$true.b2), xlab = 'True f(ROH) - bin 2', xlim = c(-0.02, 0.25), xaxt = 'n')
    # polygon(density(true.roh.bins$true.b2), col = 'grey')
    # axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    # axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    # abline(v = mean(true.roh.bins$true.b2), lty = 2, lwd = 3)
    # plot(density(true.roh.bins$true.b3), xlab = 'True f(ROH) - bin 3', xlim = c(-0.02, 0.25), xaxt = 'n')
    # polygon(density(true.roh.bins$true.b3), col = 'grey')
    # axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    # axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    # abline(v = mean(true.roh.bins$true.b3), lty = 2, lwd = 3)
    # plot(density(true.roh.bins$true.b4), xlab = 'True f(ROH) - bin 4', xlim = c(-0.02, 0.25), xaxt = 'n')
    # polygon(density(true.roh.bins$true.b4), col = 'grey')
    # axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.07)
    # axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*2, line = 0.7, lwd = 0)
    # abline(v = mean(true.roh.bins$true.b4), lty = 2, lwd = 3)
    # dev.off()

    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_hists.pdf'), width = 6, height = 6)
    figure.ct <- figure.ct + 1
    par(mar = c(6.1, 6.1, 4.1, 2.1))
    hist((true.roh.bins$true.b1), xlab = '', xlim = c(-0.02, 0.25), xaxt = 'n', ylim = c(0, 100), yaxt = 'n', ylab = '', breaks = 5, main='')
      axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.04, lwd = 3)
      axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*3, line = 2.5, lwd = 0)
      axis(2, at = c(0, 50, 100), label = rep('', 3), tck = -0.07, lwd = 3)
      axis(2, at = c(0, 100), cex.axis = text.size*3, line = 0.9, lwd = 0)
      abline(v = mean(true.roh.bins$true.b1), lty = 2, lwd = 5)
    hist((true.roh.bins$true.b2), xlab = '', xlim = c(-0.02, 0.25), xaxt = 'n', ylim = c(0, 100), yaxt = 'n', ylab = '', breaks = 5, main='')
      axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.04, lwd = 3)
      axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*3, line = 2.5, lwd = 0)
      axis(2, at = c(0, 50, 100), label = rep('', 3), tck = -0.07, lwd = 3)
      axis(2, at = c(0, 100), cex.axis = text.size*3, line = 0.9, lwd = 0)
      abline(v = mean(true.roh.bins$true.b2), lty = 2, lwd = 5)
    hist((true.roh.bins$true.b3), xlab = '', xlim = c(-0.02, 0.25), xaxt = 'n', ylim = c(0, 100), yaxt = 'n', ylab = '', breaks = 5, main='')
      axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.04, lwd = 3)
      axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*3, line = 2.5, lwd = 0)
      axis(2, at = c(0, 50, 100), label = rep('', 3), tck = -0.07, lwd = 3)
      axis(2, at = c(0, 100), cex.axis = text.size*3, line = 0.9, lwd = 0)
      abline(v = mean(true.roh.bins$true.b3), lty = 2, lwd = 5)
    hist((true.roh.bins$true.b4), xlab = '', xlim = c(-0.02, 0.25), xaxt = 'n', ylim = c(0, 100), yaxt = 'n', ylab = '', breaks = 5, main='')
      axis(1, at = c(0, 0.1, 0.2), label = rep('', 3), tck = -0.04, lwd = 3)
      axis(1, at = c(0, 0.1, 0.2), cex.axis = text.size*3, line = 2.5, lwd = 0)
      axis(2, at = c(0, 50, 100), label = rep('', 3), tck = -0.07, lwd = 3)
      axis(2, at = c(0, 100), cex.axis = text.size*3, line = 0.9, lwd = 0)
      abline(v = mean(true.roh.bins$true.b4), lty = 2, lwd = 5)
    dev.off()


    # text.size <- 1.75
    # lwd <- 2
    # pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_combined_hists.pdf'), width = 8, height = 5.5)
    # hist(true.roh.bins$true.b1, xlim = c(0, 0.2), ylim = c(0, 100), breaks = 10, cex.axis = text.size, cex.lab = text.size)
    #   hist(true.roh.bins$true.b2, add = TRUE, breaks = 10)
    #   hist(true.roh.bins$true.b3, add = TRUE, breaks = 10)
    #   hist(true.roh.bins$true.b4, add = TRUE, breaks = 10)
    #
    # hist(true.roh.bins$true.b4, xlim = c(0, 0.2), ylim = c(0, 100), breaks = 10, cex.axis = text.size, cex.lab = text.size, col = alpha('grey20', 0.5))
    #   hist(true.roh.bins$true.b3, add = TRUE, breaks = 10, col = alpha('grey40', 0.5))
    #   hist(true.roh.bins$true.b2, add = TRUE, breaks = 10, col = alpha('grey60', 0.5))
    #   hist(true.roh.bins$true.b1, add = TRUE, breaks = 10, col = alpha('grey80', 0.5))
    # dev.off()

    ##### 95. Plotting true ROH length vs. # of called ROHs overlapping #####
    alph <- 0.25
    y.max <- max(max(table(gt.overlap$called.roh.id)),
                 max(table(pl.overlap$called.roh.id)),
                 max(table(plink.overlap$called.roh.id)))
    x.max <- max(gt.overlap$len, pl.overlap$len, plink.overlap$len)
    txt.size <- 1.25
    pt.sz <- 1.25
    ylab <- 0.8

    OUT1 <- NULL
    t <- 1
    while(t == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_truelength_vs_numcalled.pdf'), width=14, height=10)
      figure.ct <- figure.ct + 1
      par(mfrow=c(3,6), mar=c(4.1,1.1,2.1,2.1))

      ## GT
      for(c in sort(unique(gt.overlap$covg))){
        print(c)
        if(c == min(gt.overlap$covg)){
          plot.new()
          # legend('right', legend='Genotypes only', bty='n', text.font=2, cex=1.25)
          mtext(side = 4, las = 3, text = 'Genotypes only:\nNumber of true ROHs overlapping\nsingle called ROH\n', cex = ylab)
        }
        OUT <- NULL
        for(i in unique(gt.overlap$id)){
          sub <- gt.overlap[gt.overlap$covg == c & gt.overlap$id == i,]
          for(r in unique(sub$called.roh.id)){ ## for each called ROH that overlaps a true ROH,
            save <- c(i, sub[sub$called.roh.id == r, 'len'][1], length(unique(sub[sub$called.roh.id == r, 'true.roh.id'])))
            save1 <- c(i, c, sub[sub$called.roh.id == r, 'len'][1], length(unique(sub[sub$called.roh.id == r, 'true.roh.id'])), 'GT')
            OUT <- rbind(OUT, save)
            OUT1 <- rbind(OUT1, save1)
          }
        }
        plot(1,1, col='transparent', xlim=c(1e5, x.max), xlab='Called ROH length (Mb)', main=paste0(c,'X'), ylim = c(1, y.max),
             cex.axis = txt.size, cex.lab = txt.size, xaxt = 'n', ylab = '')
        points(OUT[,2], OUT[,3], col = alpha(gt.col, alph), pch = 16, cex = pt.sz)
        axis(1, at = seq(from = 0, to = x.max, by = 1e6), labels = seq(from = 0, to = x.max, by = 1e6)/1e6, cex.axis = txt.size)
        ## for counts with >= 5 points, plot True ROH length mean +/- SD
        for(ct in as.numeric(names(table(OUT[,3])[which(table(OUT[,3]) >= 5)]))){
          lines(x = c(mean(OUT[OUT[,3] == ct, 2]), mean(OUT[OUT[,3] == ct, 2])), y = c(ct-0.25, ct+0.25), lwd = 4, col = gt.col)
        }
      }


      ## PL
      for(c in sort(unique(pl.overlap$covg))){
        print(c)
        if(c == min(pl.overlap$covg)){
          plot.new()
          # legend('right', legend='Genotype likelihoods', bty='n', text.font=2, cex=1.25)
          mtext(side = 4, las = 3, text = 'Genotype likelihoods:\nNumber of true ROHs overlapping\nsingle called ROH\n', cex = ylab)
        }
        OUT <- NULL
        for(i in unique(pl.overlap$id)){
          sub <- pl.overlap[pl.overlap$covg == c & pl.overlap$id == i,]
          for(r in unique(sub$called.roh.id)){ ## for each called ROH that overlaps a true ROH,
            save <- c(i, sub[sub$called.roh.id == r, 'len'][1], length(unique(sub[sub$called.roh.id == r, 'true.roh.id'])))
            save1 <- c(i, c, sub[sub$called.roh.id == r, 'len'][1], length(unique(sub[sub$called.roh.id == r, 'true.roh.id'])), 'PL')
            OUT <- rbind(OUT, save)
            OUT1 <- rbind(OUT1, save1)
          }
        }
        plot(1,1, col='transparent', xlim=c(1e5, x.max), xlab='Called ROH length (Mb)', main=paste0(c,'X'), ylim = c(1, y.max),
             cex.axis = txt.size, cex.lab = txt.size, xaxt = 'n', ylab = '')
        points(OUT[,2], OUT[,3], col = alpha(pl.col, alph), pch = 16, cex = pt.sz)
        axis(1, at = seq(from = 0, to = x.max, by = 1e6), labels = seq(from = 0, to = x.max, by = 1e6)/1e6, cex.axis = txt.size)
        ## for counts with >= 5 points, plot True ROH length mean +/- SD
        for(ct in as.numeric(names(table(OUT[,3])[which(table(OUT[,3]) >= 5)]))){
          lines(x = c(mean(OUT[OUT[,3] == ct, 2]), mean(OUT[OUT[,3] == ct, 2])), y = c(ct-0.25, ct+0.25), lwd = 4, col = pl.col)
        }
      }


      ## PLINK
      for(c in sort(unique(plink.overlap$covg))){
        print(c)
        if(c == min(plink.overlap$covg)){
          plot.new()
          # legend('right', legend='PLINK', bty='n', text.font=2, cex=1.25)
          mtext(side = 4, las = 3, text = 'PLINK:\nNumber of true ROHs overlapping\nsingle called ROH\n', cex = ylab)
        }
        OUT <- NULL
        for(i in unique(plink.overlap$id)){
          sub <- plink.overlap[plink.overlap$covg == c & plink.overlap$id == i,]
          for(r in unique(sub$called.roh.id)){ ## for each called ROH that overlaps a true ROH,
            save <- c(i, sub[sub$called.roh.id == r, 'len'][1], length(unique(sub[sub$called.roh.id == r, 'true.roh.id'])))
            save1 <- c(i, c, sub[sub$called.roh.id == r, 'len'][1], length(unique(sub[sub$called.roh.id == r, 'true.roh.id'])), 'PLINK')
            OUT <- rbind(OUT, save)
            OUT1 <- rbind(OUT1, save1)
          }
        }
        plot(1,1, col='transparent', xlim=c(1e5, x.max), xlab='Called ROH length (Mb)', main=paste0(c,'X'), ylim = c(1, y.max),
             cex.axis = txt.size, cex.lab = txt.size, xaxt = 'n', ylab = '')
        points(OUT[,2], OUT[,3], col = alpha(plink.col, alph), pch = 16, cex = pt.sz)
        axis(1, at = seq(from = 0, to = x.max, by = 1e6), labels = seq(from = 0, to = x.max, by = 1e6)/1e6, cex.axis = txt.size)
        ## for counts with >= 5 points, plot True ROH length mean +/- SD
        for(ct in as.numeric(names(table(OUT[,3])[which(table(OUT[,3]) >= 5)]))){
          lines(x = c(mean(OUT[OUT[,3] == ct, 2]), mean(OUT[OUT[,3] == ct, 2])), y = c(ct-0.25, ct+0.25), lwd = 4, col = plink.col)
        }
      }

      t <- t+1
      dev.off()
    }


    # ##### >>> 95A. Statistics for true ROH length vs. # of called ROHs overlapping (Poisson regression) #####
    # dat <- as.data.frame(OUT1)
    # dat[,2] <- gsub('x', '', dat[,2])
    # for(c in 1:4){
    #   dat[,c] <- as.numeric(dat[,c])
    # }
    # colnames(dat) <- c('id','covg','call.len','num.true','method')
    #
    # c <- 5
    # m <- 'PLINK'
    # par(mfrow = c(2,2))
    # ## model Poisson regression using glm()
    # pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_test.pdf'))
    # par(mfrow = c(2,2))
    # for(c in unique(dat$covg)){
    #   for(m in unique(dat$method)){
    #     sub <- dat[dat$covg == c & dat$method == m,]
    #     mod <- glm(num.true ~ call.len, sub, family = poisson(link = 'log'))
    #     summary(mod)
    #     shapiro.test(residuals(mod))
    #     plot(mod)
    #   }
    # }
    # dev.off()
    #
    # ## model linear regression using lm()
    # # for(c in unique(dat$covg)){
    # #   for(m in unique(dat$method)){
    #     sub <- dat[dat$covg == c & dat$method == m,]
    #     mod <- lm(sub$num.true ~ sub$call.len)
    #     summary(mod)
    #     shapiro.test(residuals(mod))
    #     plot(mod)
    # #   }
    # # }
    #
    # ## assign colors based on coverage level
    # pal <- colorRampPalette(c("orange", "purple4"))
    # cols <- pal(5)
    # col <- 1
    # for(c in sort(unique(dat$covg))){
    #   dat[dat$covg == c, 'temp.col'] <- cols[col]
    #   col <- col+1
    # }
    #
    # for(m in unique(dat$method)){
    #   plot(dat[dat$method == m, 'call.len'], dat[dat$method == m, 'num.true'], col = 'transparent', xlab = 'Called ROH length',
    #        ylab = 'Number of true ROHs overlapped', main = m)
    #   for(c in unique(dat$covg)){
    #     sub <- dat[dat$method == m & dat$covg == c,]
    #     abline(lm(sub$num.true ~ sub$call.len), col = sub$temp.col[1], lwd = 2)
    #     legend('topleft', lwd = 2, col = cols, legend = c('5X','10X','15X','30X','50X'))
    #   }
    # }

    ## tried sq-rt x-form, log x-form, and Poisson regression, nothing makes the assumption plots look less crazy

    ##### >>> 95B. What about multiple called ROHs : single true ROH? #####
    OUT <- NULL
    for(c in sort(unique(plink.overlap$covg))){
      sub.pl <- pl.overlap[pl.overlap$covg == c,]
      sub.gt <- gt.overlap[gt.overlap$covg == c,]
      sub.plink <- plink.overlap[plink.overlap$covg == c,]

      for(r in unique(sub.pl$true.roh.id)){
        save <- c(1, c, r, length(unique(sub.pl[sub.pl$true.roh.id == r, 'called.roh.id'])))
        OUT <- rbind(OUT, save)
      }

      for(r in unique(sub.gt$true.roh.id)){
        save <- c(2, c, r, length(unique(sub.gt[sub.gt$true.roh.id == r, 'called.roh.id'])))
        OUT <- rbind(OUT, save)
      }

      for(r in unique(sub.plink$true.roh.id)){
        save <- c(3, c, r, length(unique(sub.plink[sub.plink$true.roh.id == r, 'called.roh.id'])))
        OUT <- rbind(OUT, save)
      }
    }
    ## OUT[,1] == method; OUT[,2] == coverage; OUT[,3] == true ROH ID; OUT[,4] == # of called ROHs overlapping true ROH

    ##### >>> 95C. Creating lumping plots for combination with ROH length bin figures #####
    ## OUT1 (calculated above)
    ## [,1] == id; [,2] == coverage; [,3] == called ROH length; [,4] == # of true ROHs overlapping called ROH; [,5] == method
    overlaps <- as.data.frame(OUT1)
    colnames(overlaps) <- c('id','covg','call.len','num.true','method')
    overlaps$covg <- gsub('x', '', overlaps$covg)
    for(c in 1:4){
      overlaps[,c] <- as.numeric(overlaps[,c])
    }
    overlaps[overlaps$call.len < b1, 'bin'] <- 1
    overlaps[overlaps$call.len >= b1 & overlaps$call.len < b2, 'bin'] <- 2
    overlaps[overlaps$call.len >= b2 & overlaps$call.len < b3, 'bin'] <- 3
    overlaps[overlaps$call.len >= b3, 'bin'] <- 4
    overlaps[overlaps$method == 'PLINK', 'color'] <- plink.col
    overlaps[overlaps$method == 'PL', 'color'] <- pl.col
    overlaps[overlaps$method == 'GT', 'color'] <- gt.col
    overlaps[overlaps$method == 'GT', 'method'] <- 1
    overlaps[overlaps$method == 'PL', 'method'] <- 2
    overlaps[overlaps$method == 'PLINK', 'method'] <- 3
    overlaps$method <- as.numeric(overlaps$method)

    ## limit to 15X for simplicity for now
    alph <- 0.1
    width <- 0.2
    pt.size <- 3

    k <- 1
    while(k == 1){

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_lumping_plots.pdf'), width = 6, height = 6)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1, 5.1, 4.1, 4.1))

      overlaps <- overlaps[overlaps$covg == 15,]

      b <- 1
      sub <- overlaps[overlaps$bin == b,]
      plot(0, 0, col = 'transparent', main = b, xlim = c(0.75,3.25), ylim = c(0.8, max(sub$num.true)+0.2), xlab = '', ylab = '',
           xaxt = 'n', yaxt = 'n')
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.04)
      axis(1, at = c(1:3), cex.axis = text.size*3, line = 2, lwd = 0, label = c('G','L','P'))
      axis(2, at = c(1:2), label = rep('', 2), tck = -0.07)
      axis(2, at = c(1:2), cex.axis = text.size*3, line = 0.7, lwd = 0)
      for(m in unique(sub$method)){
        temp <- sub[sub$method == m,]
        points(temp$method, temp$num.true, col = alpha(temp$color[1], alph), pch = 16, cex = pt.size)
        lines(x = c(m-width, m+width), y = c(mean(temp$num.true), mean(temp$num.true)), lwd = 8, col = temp$color[1])
      }

      b <- 2
      sub <- overlaps[overlaps$bin == b,]
      plot(0, 0, col = 'transparent', main = b, xlim = c(0.75,3.25), ylim = c(0.8, max(sub$num.true)+0.2), xlab = '', ylab = '',
           xaxt = 'n', yaxt = 'n')
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.04)
      axis(1, at = c(1:3), cex.axis = text.size*3, line = 2, lwd = 0, label = c('G','L','P'))
      axis(2, at = c(1:3), label = rep('', 3), tck = -0.07)
      axis(2, at = c(1:3), cex.axis = text.size*3, line = 0.7, lwd = 0)
      for(m in unique(sub$method)){
        temp <- sub[sub$method == m,]
        points(temp$method, temp$num.true, col = alpha(temp$color[1], alph), pch = 16, cex = pt.size)
        lines(x = c(m-width, m+width), y = c(mean(temp$num.true), mean(temp$num.true)), lwd = 8, col = temp$color[1])
      }

      b <- 3
      sub <- overlaps[overlaps$bin == b,]
      plot(0, 0, col = 'transparent', main = b, xlim = c(0.75,3.25), ylim = c(0.8, max(sub$num.true)+0.2), xlab = '', ylab = '',
           xaxt = 'n', yaxt = 'n')
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.04)
      axis(1, at = c(1:3), cex.axis = text.size*3, line = 2, lwd = 0, label = c('G','L','P'))
      axis(2, at = c(1:5), label = rep('', 5), tck = -0.07)
      axis(2, at = c(1:5), cex.axis = text.size*3, line = 0.7, lwd = 0)
      for(m in unique(sub$method)){
        temp <- sub[sub$method == m,]
        points(temp$method, temp$num.true, col = alpha(temp$color[1], alph), pch = 16, cex = pt.size)
        lines(x = c(m-width, m+width), y = c(mean(temp$num.true), mean(temp$num.true)), lwd = 8, col = temp$color[1])
      }

      b <- 4
      sub <- overlaps[overlaps$bin == b,]
      plot(0, 0, col = 'transparent', main = b, xlim = c(0.75,3.25), ylim = c(0.8, max(sub$num.true)+0.2), xlab = '', ylab = '',
           xaxt = 'n', yaxt = 'n')
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.04)
      axis(1, at = c(1:3), cex.axis = text.size*3, line = 2, lwd = 0, label = c('G','L','P'))
      axis(2, at = c(1:9), label = rep('', 9), tck = -0.07)
      axis(2, at = c(2,4,6,8), cex.axis = text.size*3, line = 0.7, lwd = 0)
      for(m in unique(sub$method)){
        temp <- sub[sub$method == m,]
        points(temp$method, temp$num.true, col = alpha(temp$color[1], alph), pch = 16, cex = pt.size)
        lines(x = c(m-width, m+width), y = c(mean(temp$num.true), mean(temp$num.true)), lwd = 8, col = temp$color[1])
      }

      dev.off()
      k <- k+1
    }


    ##### >>> 95D. Creating lumping plots for stand-alone figure #####
    ## OUT1 (calculated above)
    ## [,1] == id; [,2] == coverage; [,3] == called ROH length; [,4] == # of true ROHs overlapping called ROH; [,5] == method
    overlaps <- as.data.frame(OUT1)
    colnames(overlaps) <- c('id','covg','call.len','num.true','method')
    overlaps$covg <- gsub('x', '', overlaps$covg)
    for(c in 1:4){
      overlaps[,c] <- as.numeric(overlaps[,c])
    }
    overlaps[overlaps$call.len < b1, 'bin'] <- 1
    overlaps[overlaps$call.len >= b1 & overlaps$call.len < b2, 'bin'] <- 2
    overlaps[overlaps$call.len >= b2 & overlaps$call.len < b3, 'bin'] <- 3
    overlaps[overlaps$call.len >= b3, 'bin'] <- 4
    overlaps[overlaps$method == 'PLINK', 'color'] <- plink.col
    overlaps[overlaps$method == 'PL', 'color'] <- pl.col
    overlaps[overlaps$method == 'GT', 'color'] <- gt.col
    overlaps[overlaps$method == 'GT', 'method'] <- 1
    overlaps[overlaps$method == 'PL', 'method'] <- 2
    overlaps[overlaps$method == 'PLINK', 'method'] <- 3
    overlaps$method <- as.numeric(overlaps$method)
    write.csv(overlaps, paste0('3_methods_results/',d,'_lumping_data_by_bins.csv'), row.names = FALSE)

    alph <- 0.05
    width <- 0.05
    pt.size <- 3
    mean.size <- 2
    text.size <- 1.5
    txt.exp <- 1
    lwd <- 2
    diffs <- c(-0.3, -0.15, 0, 0.15, 0.3)
    xlim.offset <- 0.4
    outline.col <- 'black'

    k <- 1
    while(k == 1){

      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_lumping_plots.pdf'), width = 6.2, height = 4.5)
      figure.ct <- figure.ct + 1
      par(mar = c(5.1, 5.1, 4.1, 4.1))

      b <- 1
      temp <- overlaps[overlaps$bin == b,]
      tot <- nrow(temp)

      plot(0, 0, col = 'transparent', main = b, xlim = c(1-xlim.offset, 3+xlim.offset), ylim = c(0.7, max(overlaps$num.true)+1.3), xlab = '', ylab = '',
           xaxt = 'n', cex.lab = text.size, cex.axis = text.size)
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.00)
      axis(1, at = c(1:3), cex.axis = text.size*txt.exp, line = 1.2, lwd = 0, label = c('BCFtools\nGenotypes','BCFtools\nLikelihoods','PLINK\n'))
      x <- 1
      for(c in sort(unique(temp$covg))){
        sub <- temp[temp$covg == c,]
        for(m in unique(sub$method)){
          sub1 <- sub[sub$method == m,]
          if(m == 1){
            cols <- gt.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 2){
            cols <- pl.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 3){
            cols <- plink.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
        }
        x <- x+1
      }
      legend('topright', pch = 16, col = plink.cols, pt.cex = pt.size, legend = c('5X','10X','15X','30X','50X'), cex = text.size, bty = 'n')
      legend('topright', pch = 16, col = pl.cols, pt.cex = pt.size, legend = c('','','','',''), cex = text.size, inset = c(0.17,0), bty = 'n')
      legend('topright', pch = 16, col = gt.cols, pt.cex = pt.size, legend = c('','','','',''), cex = text.size, inset = c(0.235,0), bty = 'n')

      b <- 2
      temp <- overlaps[overlaps$bin == b,]
      tot <- nrow(temp)

      plot(0, 0, col = 'transparent', main = b, xlim = c(1-xlim.offset, 3+xlim.offset), ylim = c(0.7, max(overlaps$num.true)+1.3), xlab = '', ylab = '',
           xaxt = 'n', cex.axis = text.size, cex.lab = text.size)
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.0)
      axis(1, at = c(1:3), cex.axis = text.size*txt.exp, line = 1.2, lwd = 0, label = c('BCFtools\nGenotypes','BCFtools\nLikelihoods','PLINK\n'))
      x <- 1
      for(c in sort(unique(temp$covg))){
        sub <- temp[temp$covg == c,]
        for(m in unique(sub$method)){
          sub1 <- sub[sub$method == m,]
          if(m == 1){
            cols <- gt.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 2){
            cols <- pl.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 3){
            cols <- plink.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
        }
        x <- x+1
      }

      b <- 3
      temp <- overlaps[overlaps$bin == b,]
      tot <- nrow(temp)

      plot(0, 0, col = 'transparent', main = b, xlim = c(1-xlim.offset, 3+xlim.offset), ylim = c(0.7, max(overlaps$num.true)+1.3), xlab = '', ylab = '',
           xaxt = 'n', cex.axis = text.size, cex.lab = text.size)
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.0)
      axis(1, at = c(1:3), cex.axis = text.size*txt.exp, line = 1.2, lwd = 0, label = c('BCFtools\nGenotypes','BCFtools\nLikelihoods','PLINK\n'))
      x <- 1
      for(c in sort(unique(temp$covg))){
        sub <- temp[temp$covg == c,]
        for(m in unique(sub$method)){
          sub1 <- sub[sub$method == m,]
          if(m == 1){
            cols <- gt.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 2){
            cols <- pl.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 3){
            cols <- plink.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
        }
        x <- x+1
      }

      b <- 4
      temp <- overlaps[overlaps$bin == b,]
      tot <- nrow(temp)

      plot(0, 0, col = 'transparent', main = b, xlim = c(1-xlim.offset, 3+xlim.offset), ylim = c(0.7, max(overlaps$num.true)+1.3), xlab = '', ylab = '',
           xaxt = 'n', cex.axis = text.size, cex.lab = text.size)
      axis(1, at = c(1:3), label = rep('', 3), tck = -0.0)
      axis(1, at = c(1:3), cex.axis = text.size*txt.exp, line = 1.2, lwd = 0, label = c('BCFtools\nGenotypes','BCFtools\nLikelihoods','PLINK\n'))
      x <- 1
      for(c in sort(unique(temp$covg))){
        sub <- temp[temp$covg == c,]
        for(m in unique(sub$method)){
          sub1 <- sub[sub$method == m,]
          if(m == 1){
            cols <- gt.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 2){
            cols <- pl.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
          if(m == 3){
            cols <- plink.cols
            points(sub1$method + diffs[x], sub1$num.true, col = alpha(cols[x], nrow(sub1)/tot), pch = 16, cex = pt.size)
            lines(x = c(m+diffs[x]-width, m+diffs[x]+width), y = c(mean(sub1$num.true), mean(sub1$num.true)), lwd = lwd, col = outline.col)
            points(m+diffs[x], mean(sub1$num.true), pch = 23, bg = cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          }
        }
        x <- x+1
      }

      dev.off()
      k <- k+1
    }

    ##### >>> 95E. Different lumping plot for SI showing info in different way #####
    # head(overlaps)
    # ymax <- max(max(table(gt.overlap$called.roh.id)),
    #              max(table(pl.overlap$called.roh.id)),
    #              max(table(plink.overlap$called.roh.id)))
    ymax <- 7
    txt.size <- 1.25
    ln.alph <- 0.8
    pt.size <- 3
    mean.size <- 1.5
    lwd <- 2
    outline.col <- 'black'
    pt.offset <- 0.015 ## changes inset of point overlays in legend (need to recal with dimension changes)

    k <- 1
    while(k == 1){
      ## Genotypes
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_SI_lumping_figure_panels.pdf'), width = 5, height = 5)
      figure.ct <- figure.ct + 1
      plot(0,0, xlim = c(1,4), ylim = c(1, ymax), col = 'transparent', xaxt = 'n', ylab = '', xlab = '', main = 'Genotypes',
           cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
        axis(1, at = c(1:4), cex.axis = txt.size, label = c('Short','Intermediate','Long','Very long'))
        axis(2, at = c(1:ymax), cex.axis = txt.size)
        x <- 5
        for(c in sort(unique(overlaps$covg), decreasing = TRUE)){
          sub <- overlaps[overlaps$covg == c & overlaps$method == 1,]
          lines(x = c(1:4), y = c(mean(sub[sub$bin == 1, 'num.true']),
                                  mean(sub[sub$bin == 2, 'num.true']),
                                  mean(sub[sub$bin == 3, 'num.true']),
                                  mean(sub[sub$bin == 4, 'num.true'])),
                col = alpha(gt.cols[x], ln.alph))
          points(x = c(1:4), y = c(mean(sub[sub$bin == 1, 'num.true']),
                                   mean(sub[sub$bin == 2, 'num.true']),
                                   mean(sub[sub$bin == 3, 'num.true']),
                                   mean(sub[sub$bin == 4, 'num.true'])),
                 pch = 23, bg = gt.cols[x], col = outline.col, cex = mean.size, lwd = lwd)
          x <- x-1
        }
        legend('topleft', pch = 23, pt.bg = gt.cols, col = gt.cols, cex = text.size, pt.cex = mean.size, lwd = 1, pt.lwd = lwd,
               legend = c('5X','10X','15X','30X','50X'), bty = 'n')
        legend('topleft', pch = 23, pt.bg = gt.cols, col = outline.col, cex = text.size, pt.cex = mean.size, pt.lwd = lwd,
               legend = c('','','','',''), bty = 'n', inset = c(pt.offset, 0))

      ## Likelihoods
      plot(0,0, xlim = c(1,4), ylim = c(1, ymax), col = 'transparent', xaxt = 'n', ylab = '', xlab = '', main = 'Likelihoods',
           cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(1:4), cex.axis = txt.size, label = c('Short','Intermediate','Long','Very long'))
      axis(2, at = c(1:ymax), cex.axis = txt.size)
      x <- 5
      for(c in sort(unique(overlaps$covg), decreasing = TRUE)){
        sub <- overlaps[overlaps$covg == c & overlaps$method == 2,]
        lines(x = c(1:4), y = c(mean(sub[sub$bin == 1, 'num.true']),
                                mean(sub[sub$bin == 2, 'num.true']),
                                mean(sub[sub$bin == 3, 'num.true']),
                                mean(sub[sub$bin == 4, 'num.true'])),
              col = alpha(pl.cols[x], ln.alph))
        points(x = c(1:4), y = c(mean(sub[sub$bin == 1, 'num.true']),
                                 mean(sub[sub$bin == 2, 'num.true']),
                                 mean(sub[sub$bin == 3, 'num.true']),
                                 mean(sub[sub$bin == 4, 'num.true'])),
               pch = 23, bg = pl.cols[x], col = outline.col, cex = mean.size, lwd = lwd)
        x <- x-1
      }
      legend('topleft', pch = 23, pt.bg = pl.cols, col = pl.cols, cex = text.size, pt.cex = mean.size, lwd = 1, pt.lwd = lwd,
             legend = c('5X','10X','15X','30X','50X'), bty = 'n')
      legend('topleft', pch = 23, pt.bg = pl.cols, col = outline.col, cex = text.size, pt.cex = mean.size, pt.lwd = lwd,
             legend = c('','','','',''), bty = 'n', inset = c(pt.offset, 0))

      ## PLINK
      plot(0,0, xlim = c(1,4), ylim = c(1, ymax), col = 'transparent', xaxt = 'n', ylab = '', xlab = '', main = 'PLINK',
           cex.axis = txt.size, cex.lab = txt.size, yaxt = 'n')
      axis(1, at = c(1:4), cex.axis = txt.size, label = c('Short','Intermediate','Long','Very long'))
      axis(2, at = c(1:ymax), cex.axis = txt.size)
      x <- 5
      for(c in sort(unique(overlaps$covg), decreasing = TRUE)){
        sub <- overlaps[overlaps$covg == c & overlaps$method == 3,]
        lines(x = c(1:4), y = c(mean(sub[sub$bin == 1, 'num.true']),
                                mean(sub[sub$bin == 2, 'num.true']),
                                mean(sub[sub$bin == 3, 'num.true']),
                                mean(sub[sub$bin == 4, 'num.true'])),
              col = alpha(plink.cols[x], ln.alph))
        points(x = c(1:4), y = c(mean(sub[sub$bin == 1, 'num.true']),
                                 mean(sub[sub$bin == 2, 'num.true']),
                                 mean(sub[sub$bin == 3, 'num.true']),
                                 mean(sub[sub$bin == 4, 'num.true'])),
               pch = 23, bg = plink.cols[x], col = outline.col, cex = mean.size, lwd = lwd)
        x <- x-1
      }
      legend('topleft', pch = 23, pt.bg = plink.cols, col = plink.cols, cex = text.size, pt.cex = mean.size, lwd = 1, pt.lwd = lwd,
             legend = c('5X','10X','15X','30X','50X'), bty = 'n')
      legend('topleft', pch = 23, pt.bg = plink.cols, col = outline.col, cex = text.size, pt.cex = mean.size, pt.lwd = lwd,
             legend = c('','','','',''), bty = 'n', inset = c(pt.offset, 0))

      dev.off()
      k <- k+1
    }


    ##### 96. Plotting true vs. called ROH lengths (scatterplot) #####
    alph <- 0.25
    y.max <- max(pl.overlap$len, gt.overlap$len, plink.overlap$len, true.rohs$length)

    t <- 1
    while(t == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_ROHspecificlengthcomps_scatterplot.pdf'), width=16, height=12)
      figure.ct <- figure.ct + 1
      par(mfrow=c(3,6), mar=c(4.1,3.1,2.1,2.1))
      ## PL
      for(c in sort(unique(pl.overlap$covg))){
        if(c == min(pl.overlap$covg)){
          plot.new()
          legend('right', legend='Genotype likelihoods', bty='n', text.font=2, cex=1.25)
        }
        sub <- pl.overlap[pl.overlap$covg == c,]
        plot(1,1, col='transparent', xlim=c(1e5, y.max), xlab='True length', main=paste0(c), ylim = c(1e5, y.max))
        for(i in 1:nrow(sub)){ ## for each called ROH that overlaps a true ROH, draw a line connecting its length to the length of the true ROH it overlaps
          points(true.rohs[true.rohs$true.roh.id == sub$true.roh.id[i], 'length'], sub$len[i], col=alpha(pl.col, alph), pch = 16, cex = 0.5)
        }
        abline(0,1, lty = 2, lwd = 0.5)
      }
      ## GT
      for(c in sort(unique(gt.overlap$covg))){
        if(c == min(gt.overlap$covg)){
          plot.new()
          legend('right', legend='Genotypes only', bty='n', text.font=2, cex=1.25)
        }
        sub <- gt.overlap[gt.overlap$covg == c,]
        plot(1,1, col='transparent', xlim=c(1e5, y.max), xlab='True length', main=paste0(c), ylim = c(1e5, y.max))
        for(i in 1:nrow(sub)){ ## for each called ROH that overlaps a true ROH, draw a line connecting its length to the length of the true ROH it overlaps
          points(true.rohs[true.rohs$true.roh.id == sub$true.roh.id[i], 'length'], sub$len[i], col=alpha(gt.col, alph), pch = 16, cex = 0.5)
        }
        abline(0,1, lty = 2, lwd = 0.5)
      }
      ## PLINK
      for(c in sort(unique(plink.overlap$covg))){
        if(c == min(plink.overlap$covg)){
          plot.new()
          legend('right', legend='PLINK', bty='n', text.font=2, cex=1.25)
        }
        sub <- plink.overlap[plink.overlap$covg == c,]
        plot(1,1, col='transparent', xlim=c(1e5, y.max), xlab='True length', main=paste0(c), ylim = c(1e5, y.max))
        for(i in 1:nrow(sub)){ ## for each called ROH that overlaps a true ROH, draw a line connecting its length to the length of the true ROH it overlaps
          points(true.rohs[true.rohs$true.roh.id == sub$true.roh.id[i], 'length'], sub$len[i], col=alpha(plink.col, alph), pch = 16, cex = 0.5)
        }
        abline(0,1, lty = 2, lwd = 0.5)
      }
      dev.off()
      t <- t+1
    }

    ##### 97. Plotting true vs. called ROH lengths (line plot - shows lumping pattern in calls) #####
    alph <- 0.2
    y.max <- max(pl.overlap$len, gt.overlap$len, plink.overlap$len, true.rohs$length)

    t <- 1
    while(t == 1){
      pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_ROHspecificlengthcomps_lineplot.pdf'), width=16, height=12)
      figure.ct <- figure.ct + 1
      par(mfrow=c(3,6), mar=c(4.1,3.1,2.1,2.1))
      ## PL
      for(c in sort(unique(pl.overlap$covg))){
        if(c == min(pl.overlap$covg)){
          plot.new()
          legend('right', legend='Genotype likelihoods', bty='n', text.font=2, cex=1.25)
        }
        sub <- pl.overlap[pl.overlap$covg == c,]
        plot(1,1, col='transparent', xlim=c(1,2), xaxt = 'n', xlab='', main=paste0(c), ylim = c(1e5, y.max))
        for(i in 1:nrow(sub)){ ## for each called ROH that overlaps a true ROH, draw a line connecting its length to the length of the true ROH it overlaps
          lines(c(1,2), c(true.rohs[true.rohs$true.roh.id == sub$true.roh.id[i], 'length'], sub$len[i]), col=alpha(pl.col, alph))
        }
        axis(1, at=c(1,2), labels=c('True','Called'), padj = 0.5)
      }
      ## GT
      for(c in sort(unique(gt.overlap$covg))){
        if(c == min(gt.overlap$covg)){
          plot.new()
          legend('right', legend='Genotypes only', bty='n', text.font=2, cex=1.25)
        }
        sub <- gt.overlap[gt.overlap$covg == c,]
        plot(1,1, col='transparent', xlim=c(1,2), xaxt = 'n', xlab='', main=paste0(c), ylim = c(1e5, y.max))
        for(i in 1:nrow(sub)){ ## for each called ROH that overlaps a true ROH, draw a line connecting its length to the length of the true ROH it overlaps
          lines(c(1,2), c(true.rohs[true.rohs$true.roh.id == sub$true.roh.id[i], 'length'], sub$len[i]), col=alpha(gt.col, alph))
        }
        axis(1, at=c(1,2), labels=c('True','Called'), padj = 0.5)
      }
      ## PLINK
      for(c in sort(unique(plink.overlap$covg))){
        if(c == min(plink.overlap$covg)){
          plot.new()
          legend('right', legend='PLINK', bty='n', text.font=2, cex=1.25)
        }
        sub <- plink.overlap[plink.overlap$covg == c,]
        plot(1,1, col='transparent', xlim=c(1,2), xaxt = 'n', xlab='', main=paste0(c), ylim = c(1e5, y.max))
        for(i in 1:nrow(sub)){ ## for each called ROH that overlaps a true ROH, draw a line connecting its length to the length of the true ROH it overlaps
          lines(c(1,2), c(true.rohs[true.rohs$true.roh.id == sub$true.roh.id[i], 'length'], sub$len[i]), col=alpha(plink.col, alph))
        }
        axis(1, at=c(1,2), labels=c('True','Called'), padj = 0.5)
      }
      dev.off()
      t <- t+1
    }


    ##### 98. Plotting false-positives and -negatives as a function of true f(ROH) #####
    pal <- colorRampPalette(c("orange", "purple4"))
    cols <- pal(5)
    # cols <- brewer.pal(5, 'Spectral')

    plink.true.v.called[plink.true.v.called$covg == 5, 'temp.colour'] <- cols[1]
    gt.true.v.called[gt.true.v.called$covg == 5, 'temp.colour'] <- cols[1]
    pl.true.v.called[pl.true.v.called$covg == 5, 'temp.colour'] <- cols[1]
    plink.true.v.called[plink.true.v.called$covg == 10, 'temp.colour'] <- cols[2]
    gt.true.v.called[gt.true.v.called$covg == 10, 'temp.colour'] <- cols[2]
    pl.true.v.called[pl.true.v.called$covg == 10, 'temp.colour'] <- cols[2]
    plink.true.v.called[plink.true.v.called$covg == 15, 'temp.colour'] <- cols[3]
    gt.true.v.called[gt.true.v.called$covg == 15, 'temp.colour'] <- cols[3]
    pl.true.v.called[pl.true.v.called$covg == 15, 'temp.colour'] <- cols[3]
    plink.true.v.called[plink.true.v.called$covg == 30, 'temp.colour'] <- cols[4]
    gt.true.v.called[gt.true.v.called$covg == 30, 'temp.colour'] <- cols[4]
    pl.true.v.called[pl.true.v.called$covg == 30, 'temp.colour'] <- cols[4]
    plink.true.v.called[plink.true.v.called$covg == 50, 'temp.colour'] <- cols[5]
    gt.true.v.called[gt.true.v.called$covg == 50, 'temp.colour'] <- cols[5]
    pl.true.v.called[pl.true.v.called$covg == 50, 'temp.colour'] <- cols[5]

    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_true_fROH_vs_false_pos_and_neg_rates.pdf'), width = 7, height = 6)
    figure.ct <- figure.ct + 1
    par(xpd = TRUE, mar = c(5.1, 4.1, 4.1, 4.1))
    plot(plink.true.v.called$true.froh, plink.true.v.called$false.pos.rate, pch = 16, col = plink.true.v.called$temp.colour, xlab = 'True f(ROH)', ylab = 'False positive rate', main = 'PLINK')
    legend('topright', col = cols, legend = c('5X','10X','15X','30X','50X'), pch = 16, inset = c(-0.15, 0))
    plot(pl.true.v.called$true.froh, pl.true.v.called$false.pos.rate, pch = 16, col = pl.true.v.called$temp.colour, xlab = 'True f(ROH)', ylab = 'False positive rate', main = 'Genotype likelihoods')
    legend('topright', col = cols, legend = c('5X','10X','15X','30X','50X'), pch = 16, inset = c(-0.15, 0))
    plot(gt.true.v.called$true.froh, gt.true.v.called$false.pos.rate, pch = 16, col = gt.true.v.called$temp.colour, xlab = 'True f(ROH)', ylab = 'False positive rate', main = 'Genotypes only')
    legend('topright', col = cols, legend = c('5X','10X','15X','30X','50X'), pch = 16, inset = c(-0.15, 0))


    plot(plink.true.v.called$true.froh, plink.true.v.called$false.neg.rate, pch = 16, col = plink.true.v.called$temp.colour, xlab = 'True f(ROH)', ylab = 'False negative rate', main = 'PLINK')
    legend('topright', col = cols, legend = c('5X','10X','15X','30X','50X'), pch = 16, inset = c(-0.15, 0))
    plot(pl.true.v.called$true.froh, pl.true.v.called$false.neg.rate, pch = 16, col = gt.true.v.called$temp.colour, xlab = 'True f(ROH)', ylab = 'False negative rate', main = 'Genotype likelihoods')
    legend('topright', col = cols, legend = c('5X','10X','15X','30X','50X'), pch = 16, inset = c(-0.15, 0))
    plot(gt.true.v.called$true.froh, gt.true.v.called$false.neg.rate, pch = 16, col = pl.true.v.called$temp.colour, xlab = 'True f(ROH)', ylab = 'False negative rate', main = 'Genotypes only')
    legend('topright', col = cols, legend = c('5X','10X','15X','30X','50X'), pch = 16, inset = c(-0.15, 0))
    dev.off()

    ##### 99. Plotting individual true ROHs and called ROHs for all 3 analyses and coverage levels #####
    wid <- 3 ## line widths in plot ~ or ~
    hit <- 0.1 ## polygon height
    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_individual_true_and_called_ROHs_across_chromosome.pdf'), width = 15, height = 7)
    figure.ct <- figure.ct + 1
    par(mar = c(5.1, 5.1, 4.1, 2.1))
    for(i in sort(as.numeric(unique(bcf.pl.res$id)))){
      sub.pl <- bcf.pl.res[bcf.pl.res$id == i,]
      sub.gt <- bcf.gt.res[bcf.gt.res$id == i,]
      sub.pk <- plink.res[plink.res$id == i,]
      sub.true <- true.rohs[true.rohs$id == i,]
      plot(0,0, xlim = c(1, chrom.len), ylim = c(1, 16), col = 'transparent', yaxt = 'n', xlab = 'Chromosome position (bp)', ylab = '', main = i)
      lines(x = c(-2e6, 40e6), y = c(1.5, 1.5), lty = 2)
      lines(x = c(-2e6, 40e6), y = c(6.5, 6.5), lty = 2)
      lines(x = c(-2e6, 40e6), y = c(11.5, 11.5), lty = 2)
      axis(2, at = c(1:16), labels = c('True','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X'), las = 2)
      for(r in 1:nrow(sub.true)){
        # lines(x = c(sub.true$start[r], sub.true$end[r]), y = c(1, 1), col = 'black', lwd = wid)
        polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1, 1, 18, 18),
                col = alpha('lightgrey', 0.3), border = NA)
        polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1-hit, 1-hit, 1+hit, 1+hit),
                col = 'black', border = NA)
      }
      y <- 2
      for(c in c(5, 10, 15, 30, 50)){
        temp <- sub.pl[sub.pl$covg == c,]
        if(nrow(temp) > 0){
          for(r in 1:nrow(temp)){
            # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = pl.col, lwd = wid)
            polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                    col = pl.col, border = NA)
          }
        }
        y <- y+1
      }
      for(c in c(5, 10, 15, 30, 50)){
        temp <- sub.gt[sub.gt$covg == c,]
        if(nrow(temp) > 0){
          for(r in 1:nrow(temp)){
            # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = gt.col, lwd = wid)
            polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                    col = gt.col, border = NA)
          }
        }
        y <- y+1
      }
      for(c in c(5, 10, 15, 30, 50)){
        temp <- sub.pk[sub.pk$covg == c,]
        if(nrow(temp) > 0){
          for(r in 1:nrow(temp)){
            # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = plink.col, lwd = wid)
            polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                    col = plink.col, border = NA)
          }
        }
        y <- y+1
      }
      mtext('PLINK', side = 2, line = 3, at = 14)
      mtext('BCFtools\nGenotypes', side = 2, line = 3, at = 9)
      mtext('BCFtools\nLikelihoods', side = 2, line = 3, at = 4)
    }
    dev.off()

    ##### >>> 99A. Plotting individual true ROHs and called ROHs for all 3 analyses and coverage levels - specific individual and window #####
    wid <- 3 ## line widths in plot ~ or ~
    hit <- 0.1 ## polygon height
    xmin <- 1.25e7
    xmax <- 2.0e7

    pdf(paste0('/Users/Avril/Desktop/',d,'_read_dat_sample.pdf'), width = 15, height = 7)
    par(mar = c(5.1, 5.1, 4.1, 2.1))
    for(i in 26){
      sub.pl <- bcf.pl.res[bcf.pl.res$id == i,]
      sub.gt <- bcf.gt.res[bcf.gt.res$id == i,]
      sub.pk <- plink.res[plink.res$id == i,]
      sub.true <- true.rohs[true.rohs$id == i,]
      plot(0,0, ylim = c(1, 16), col = 'transparent', yaxt = 'n', xlab = 'Chromosome position (bp)', ylab = '', main = i,
           xlim = c(xmin, xmax))
      lines(x = c(-2e6, 40e6), y = c(1.5, 1.5), lty = 2)
      lines(x = c(-2e6, 40e6), y = c(6.5, 6.5), lty = 2)
      lines(x = c(-2e6, 40e6), y = c(11.5, 11.5), lty = 2)
      axis(2, at = c(1:16), labels = c('True','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X'), las = 2)
      for(r in 1:nrow(sub.true)){
        # lines(x = c(sub.true$start[r], sub.true$end[r]), y = c(1, 1), col = 'black', lwd = wid)
        polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1, 1, 18, 18),
                col = alpha('lightgrey', 0.3), border = NA)
        polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1-hit, 1-hit, 1+hit, 1+hit),
                col = 'black', border = NA)
      }
      y <- 2
      for(c in c(5, 10, 15, 30, 50)){
        temp <- sub.pl[sub.pl$covg == c,]
        if(nrow(temp) > 0){
          for(r in 1:nrow(temp)){
            # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = pl.col, lwd = wid)
            polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                    col = pl.col, border = NA)
          }
        }
        y <- y+1
      }
      for(c in c(5, 10, 15, 30, 50)){
        temp <- sub.gt[sub.gt$covg == c,]
        if(nrow(temp) > 0){
          for(r in 1:nrow(temp)){
            # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = gt.col, lwd = wid)
            polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                    col = gt.col, border = NA)
          }
        }
        y <- y+1
      }
      for(c in c(5, 10, 15, 30, 50)){
        temp <- sub.pk[sub.pk$covg == c,]
        if(nrow(temp) > 0){
          for(r in 1:nrow(temp)){
            # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = plink.col, lwd = wid)
            polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                    col = plink.col, border = NA)
          }
        }
        y <- y+1
      }
      mtext('PLINK', side = 2, line = 3, at = 14)
      mtext('BCFtools\nGenotypes', side = 2, line = 3, at = 9)
      mtext('BCFtools\nLikelihoods', side = 2, line = 3, at = 4)
    }
    dev.off()
  }

  ##### Plots comparing default and viterbi results for BCFtools #####
  figure.ct <- 37
  ##### 1. Read in data and summary files #####
  ### True ROH information
  true.rohs <- read.table(paste0('slim_true_data/',d,'_true_roh_coords.txt'))
  colnames(true.rohs) <- c('id','start','end','length','true.roh.id')
  true.rohs <- true.rohs[true.rohs$length >= 100000,]
  chrom.len <- 30e6
  if(d == 'decline'){
    true.rohs <- true.rohs[true.rohs$id %notin% c(33,48),]
  }

  ### Read in ROH calling results
  bcf.res <- read.table('bcftools_output/bcftoolsROH_all_coordinates.txt',
                        header = TRUE, sep = '\t')
  bcf.res$covg <- as.numeric(gsub('x', '', bcf.res$covg))
  bcf.res$id <- gsub('i', '', bcf.res$id)
  bcf.res <- bcf.res[bcf.res$demo == d,]
  if(d == 'decline'){
    bcf.res <- bcf.res[bcf.res$id %notin% c(33,48),]
  }
  bcf.res <- bcf.res[bcf.res$length >= 100000,] ## applying 100kb filter

  ## GT
  bcf.gt.def.res <- bcf.res[bcf.res$method == 'GT' & bcf.res$hmm == 'defaults',]
  bcf.gt.vit.res <- bcf.res[bcf.res$method == 'GT' & bcf.res$hmm == 'vtrained',]

  ## PL
  bcf.pl.def.res <- bcf.res[bcf.res$method == 'PL' & bcf.res$hmm == 'defaults',]
  bcf.pl.vit.res <- bcf.res[bcf.res$method == 'PL' & bcf.res$hmm == 'vtrained',]

  ### Write/read in called vs. true f(ROH) results
  OUT <- NULL
  for(i in unique(true.rohs$id)){
    true.froh <- sum(true.rohs[true.rohs$id == i, 'length'])/chrom.len
    for(c in unique(bcf.gt.vit.res$covg)){
      gt.def.froh <- sum(bcf.gt.def.res[bcf.gt.def.res$id == i & bcf.gt.def.res$covg == c, 'length'])/chrom.len
      gt.vit.froh <- sum(bcf.gt.vit.res[bcf.gt.vit.res$id == i & bcf.gt.vit.res$covg == c, 'length'])/chrom.len
      pl.def.froh <- sum(bcf.pl.def.res[bcf.pl.def.res$id == i & bcf.pl.def.res$covg == c, 'length'])/chrom.len
      pl.vit.froh <- sum(bcf.pl.vit.res[bcf.pl.vit.res$id == i & bcf.pl.vit.res$covg == c, 'length'])/chrom.len
      save <- c(i, c, true.froh, gt.def.froh, gt.vit.froh, pl.def.froh, pl.vit.froh)
      OUT <- rbind(OUT, save)
    }
  }
  colnames(OUT) <- c('id','covg','true.froh','gt.def.froh', 'gt.vit.froh', 'pl.def.froh', 'pl.vit.froh')
  write.csv(OUT, paste0('bcftools_output/',d,'_true_vs_called_froh_data.csv'), row.names = FALSE)
  froh.stats <- read.csv(paste0('bcftools_output/',d,'_true_vs_called_froh_data.csv'), header = TRUE)

  ### Read in false pos/neg rate information
  pl.def.true.v.called <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_defaults_PL.csv'))
  pl.vit.true.v.called <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_vtrained_PL.csv'))
  gt.def.true.v.called <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_defaults_GT.csv'))
  gt.vit.true.v.called <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_vtrained_GT.csv'))
  
  print(paste0(d,' - ',max(pl.def.true.v.called$false.pos.rate, pl.vit.true.v.called$false.pos.rate,
                           gt.def.true.v.called$false.pos.rate, gt.vit.true.v.called$false.pos.rate)))

  ##### 2. Plot false positive and negative rates #####
  i <- 1
  while(i == 1){ ## a loop just to make all the separate plots write with 1 command
    text.size <- 1.25
    shrink <- 2000 ## higher # here ==> narrower x-direction spread for points
    alph <- 0.4
    poly.alph <- 0.5
    pt.cex <- 0.7
    false.pos.y.lim <- 0.1
    percentile <- 50
    low.q <- (100-percentile)/100/2
    hi.q <- 1 - (100-percentile)/100/2
    diff <- 0.25 ## distance between defaults and vtrained
    OUT <- NULL ## for storing summary statistics, false pos
    OUT1 <- NULL ## false neg
    
    pl.pal <- colorRampPalette(c(lt.pl.col, max.covg.col))
    pl.cols <- pl.pal(5)
    gt.pal <- colorRampPalette(c(gt.col, max.covg.col))
    gt.cols <- gt.pal(5)

    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_bcftools_viterbi_v_defaults_false_pos_neg_rates.pdf'), width = 9, height = 12)
    figure.ct <- figure.ct + 1
    par(mar = c(5.1,6.1,6.1,2.1), mfrow = c(2,1))
    plot(pl.def.true.v.called$covg, pl.def.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim),
         xlab = 'Coverage', ylab = 'False positive rate', col = 'transparent',
         xlim = c(0.65,10.35), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
    mtext(scen.name, side = 3, adj = 0, font = 3, cex = text.size, line = 0.25)
    axis(2, at = seq(0, false.pos.y.lim, 0.025), labels = c('0.00','','0.05','','0.10'), cex.axis = text.size)
    axis(1, at = seq(1.5, 10.5, 2), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
    x <- 0.75
    for(c in unique(gt.def.true.v.called$covg)){
      temp.def <- gt.def.true.v.called[gt.def.true.v.called$covg == c,]
      temp.vit <- gt.vit.true.v.called[gt.vit.true.v.called$covg == c,]

      ## defaults
      f <- sample(c(-100:100), length(unique(temp.def$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(gt.cols[3], poly.alph), border = NA)

      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.def$false.pos.rate, na.rm = TRUE), median(temp.def$false.pos.rate, na.rm = TRUE)),
            pch = 16, col = gt.cols[3], lwd = 2)
      points(f, temp.def$false.pos.rate, pch = 16, col = alpha(gt.cols[3], alph), cex = pt.cex)
      OUT <- rbind(OUT, c(d,'GT','defaults', median(temp.def$false.pos.rate), mean(temp.def$false.pos.rate), sd(temp.def$false.pos.rate)))
      x <- x+0.5
      
      ## vtrained
      f <- sample(c(-100:100), length(unique(temp.vit$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(gt.cols[1], poly.alph), border = NA)

      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.vit$false.pos.rate, na.rm = TRUE), median(temp.vit$false.pos.rate, na.rm = TRUE)),
            pch = 16, col = gt.cols[1], lwd = 2)
      points(f, temp.vit$false.pos.rate, pch = 16, col = alpha(gt.cols[1], alph), cex = pt.cex)
      OUT <- rbind(OUT, c(d,'GT','vtrained', median(temp.vit$false.pos.rate), mean(temp.vit$false.pos.rate), sd(temp.vit$false.pos.rate)))
      x <- x+1.5
    }

    x <- 1.75
    for(c in unique(pl.def.true.v.called$covg)){
      temp.def <- pl.def.true.v.called[pl.def.true.v.called$covg == c,]
      temp.vit <- pl.vit.true.v.called[pl.vit.true.v.called$covg == c,]
      
      ## defaults
      f <- sample(c(-100:100), length(unique(temp.def$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.def$false.pos.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(pl.cols[3], poly.alph), border = NA)

      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.def$false.pos.rate, na.rm = TRUE), median(temp.def$false.pos.rate, na.rm = TRUE)),
            pch = 16, col = pl.cols[3], lwd = 2)
      points(f, temp.def$false.pos.rate, pch = 16, col = alpha(pl.cols[3], alph), cex = pt.cex)
      OUT <- rbind(OUT, c(d,'PL','defaults', median(temp.def$false.pos.rate), mean(temp.def$false.pos.rate), sd(temp.def$false.pos.rate)))
      x <- x+0.5

      ## vtrained
      f <- sample(c(-100:100), length(unique(temp.vit$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.vit$false.pos.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(pl.cols[1], poly.alph), border = NA)

      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.vit$false.pos.rate, na.rm = TRUE), median(temp.vit$false.pos.rate, na.rm = TRUE)),
            pch = 16, col = pl.cols[1], lwd = 2)
      points(f, temp.vit$false.pos.rate, pch = 16, col = alpha(pl.cols[1], alph), cex = pt.cex)
      OUT <- rbind(OUT, c(d,'PL','vtrained', median(temp.vit$false.pos.rate), mean(temp.vit$false.pos.rate), sd(temp.vit$false.pos.rate)))
      x <- x+1.5
    }

    ## turn on or off vertical dashed lines between coverage levels
    abline(v = 2.5, lty = 2, col = 'grey')
    abline(v = 4.5, lty = 2, col = 'grey')
    abline(v = 6.5, lty = 2, col = 'grey')
    abline(v = 8.5, lty = 2, col = 'grey')
    
    par(xpd = TRUE)
    legend('top', fill = c(gt.cols[3], gt.cols[1], pl.cols[3], pl.cols[1]), ncol = 2,
           legend = c('Genotypes - Defaults','Genotypes - Viterbi-trained','Likelihoods - Defaults','Likelihoods - Viterbi-trained'),
           inset = c(0, -0.3), bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)

    par(xpd = FALSE)
    plot(pl.def.true.v.called$covg, pl.def.true.v.called$false.neg.rate, ylim = c(0, 1),
         xlab = 'Coverage', ylab = 'False negative rate', col = 'transparent',
         xlim = c(0.65,10.35), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
    axis(2, at = seq(0, 1, 0.25), cex.axis = text.size)
    axis(1, at = seq(1.5, 10.5, 2), labels = c('5X','10X','15X','30X','50X'), cex.axis = text.size)
    x <- 0.75
    for(c in unique(gt.def.true.v.called$covg)){
      temp.def <- gt.def.true.v.called[gt.def.true.v.called$covg == c,]
      temp.vit <- gt.vit.true.v.called[gt.vit.true.v.called$covg == c,]
      
      ## defaults
      f <- sample(c(-100:100), length(unique(temp.def$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(gt.cols[3], poly.alph), border = NA)
      
      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.def$false.neg.rate, na.rm = TRUE), median(temp.def$false.neg.rate, na.rm = TRUE)),
            pch = 16, col = gt.cols[3], lwd = 2)
      points(f, temp.def$false.neg.rate, pch = 16, col = alpha(gt.cols[3], alph), cex = pt.cex)
      OUT1 <- rbind(OUT1, c(d,'GT','defaults', median(temp.def$false.neg.rate), mean(temp.def$false.neg.rate), sd(temp.def$false.neg.rate)))
      x <- x+0.5
      
      ## vtrained
      f <- sample(c(-100:100), length(unique(temp.vit$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(gt.cols[1], poly.alph), border = NA)
      
      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.vit$false.neg.rate, na.rm = TRUE), median(temp.vit$false.neg.rate, na.rm = TRUE)),
            pch = 16, col = gt.cols[1], lwd = 2)
      points(f, temp.vit$false.neg.rate, pch = 16, col = alpha(gt.cols[1], alph), cex = pt.cex)
      OUT1 <- rbind(OUT1, c(d,'GT','vtrained', median(temp.vit$false.neg.rate), mean(temp.vit$false.neg.rate), sd(temp.vit$false.neg.rate)))
      x <- x+1.5
    }
    
    x <- 1.75
    for(c in unique(pl.def.true.v.called$covg)){
      temp.def <- pl.def.true.v.called[pl.def.true.v.called$covg == c,]
      temp.vit <- pl.vit.true.v.called[pl.vit.true.v.called$covg == c,]
      
      ## defaults
      f <- sample(c(-100:100), length(unique(temp.def$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.def$false.neg.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(pl.cols[3], poly.alph), border = NA)
      
      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.def$false.neg.rate, na.rm = TRUE), median(temp.def$false.neg.rate, na.rm = TRUE)),
            pch = 16, col = pl.cols[3], lwd = 2)
      points(f, temp.def$false.neg.rate, pch = 16, col = alpha(pl.cols[3], alph), cex = pt.cex)
      OUT1 <- rbind(OUT1, c(d,'PL','defaults', median(temp.def$false.neg.rate), mean(temp.def$false.neg.rate), sd(temp.def$false.neg.rate)))
      x <- x+0.5
      
      ## vtrained
      f <- sample(c(-100:100), length(unique(temp.vit$id)))
      f <- f/shrink+x
      polygon(x = c(x-diff/2, x+diff/2, x+diff/2, x-diff/2), y = c(quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp.vit$false.neg.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(pl.cols[1], poly.alph), border = NA)
      
      lines(x = c(x-diff/2, x+diff/2), y = c(median(temp.vit$false.neg.rate, na.rm = TRUE), median(temp.vit$false.neg.rate, na.rm = TRUE)),
            pch = 16, col = pl.cols[1], lwd = 2)
      points(f, temp.vit$false.neg.rate, pch = 16, col = alpha(pl.cols[1], alph), cex = pt.cex)
      OUT1 <- rbind(OUT1, c(d,'PL','vtrained', median(temp.vit$false.neg.rate), mean(temp.vit$false.neg.rate), sd(temp.vit$false.neg.rate)))
      x <- x+1.5
    }
    
    ## turn on or off vertical dashed lines between coverage levels
    abline(v = 2.5, lty = 2, col = 'grey')
    abline(v = 4.5, lty = 2, col = 'grey')
    abline(v = 6.5, lty = 2, col = 'grey')
    abline(v = 8.5, lty = 2, col = 'grey')
    
    # par(xpd = TRUE)
    # legend('top', fill = c(gt.col, lt.gt.col, pl.col, lt.pl.col), ncol = 2,
    #        legend = c('Genotypes - Defaults','Genotypes - Viterbi-trained','Likelihoods - Defaults','Likelihoods - Viterbi-trained'),
    #        inset = c(0, -0.25), bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)
    
    dev.off()
    
    colnames(OUT) <- c('demo','method','hmm','median','mean','sd')
    write.csv(OUT, paste0('bcftools_output/',d,'_viterbi_v_default_false_pos_summary_statistics.csv'),
              row.names = FALSE, quote = FALSE)
    colnames(OUT1) <- c('demo','method','hmm','median','mean','sd')
    write.csv(OUT1, paste0('bcftools_output/',d,'_viterbi_v_default_false_neg_summary_statistics.csv'),
              row.names = FALSE, quote = FALSE)
    
    i <- i+1
  }
  ##### 3. Calculate some stats for comparing true vs. called f(ROH) relationships across methods #####
  sets <- cbind(c('GT_defaults','GT_vtrained','PL_defaults','PL_vtrained'), c(4:7))
  for(r in 1:nrow(sets)){
    colm <- as.numeric(sets[r,2])
    name <- sets[r,1]
    
    OUT <- NULL ## for saving statsy info
    
    pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_','viterbi_vs_default_truefROH_vs_callfROH_lm_assumptions_',name,'.pdf'), 
        width = 10, height = 10)
    figure.ct <- figure.ct + 1
    par(mfrow = c(2,2))
    
    for(c in sort(unique(froh.stats$covg))){
      print(paste0(name,' - ',c))
      sub <- froh.stats[froh.stats$covg == c,]
      mod <- lm(sub[,colm] ~ sub$true.froh)
      # mod <- lm(log(sub$call.froh + 1) ~ sub$true.froh) ## trying log x-form to deal with heteroscedasticity in PL models (didn't help)
      int.lo.lim <- confint(mod, level = 0.95)[1,1]
      int.up.lim <- confint(mod, level = 0.95)[1,2]
      slope.lo.lim <- confint(mod, level = 0.95)[2,1]
      slope.up.lim <- confint(mod, level = 0.95)[2,2]
      int <- mod$coefficients[1]
      slope <- mod$coefficients[2]
      save <- c(c, summary(mod)$adj.r.squared, summary(mod)$coefficients[2,4], 
                int.lo.lim, int, int.up.lim, slope.lo.lim, slope, slope.up.lim)
      OUT <- rbind(OUT, save)
  
      ## test for normality of residuals (none print, so all good here)
      if(shapiro.test(residuals(mod))[2] < 0.05){
        print(paste0('Residuals not normal: ',name,' - ',c,'X'))
      }
      plot(mod, main = paste0(name,' - ',c,'X'))
  
      ## test whether intercept differs from 0 using 95% CIs
      if(confint(mod, level = 0.95)[1,1] <= 0 & confint(mod, level = 0.95)[1,2] >= 0){
        # print(paste0(name,' - ',c,'X - Intercept == 0'))
      }
      if((confint(mod, level = 0.95)[1,1] <= 0 & confint(mod, level = 0.95)[1,2] <= 0) |
         (confint(mod, level = 0.95)[1,1] >= 0 & confint(mod, level = 0.95)[1,2] >= 0)){
        # print(paste0(name,' - ',c,'X - Intercept != 0'))
      }
  
      ## test whether slopes differ from 1 using 95% CIs
      if(confint(mod, level = 0.95)[2,1] <= 1 & confint(mod, level = 0.95)[2,2] >= 1){
        # print(paste0(name,' - ',c,'X - Slope == 1'))
      }
      if((confint(mod, level = 0.95)[2,1] <= 1 & confint(mod, level = 0.95)[2,2] <= 1) |
         (confint(mod, level = 0.95)[2,1] >= 1 & confint(mod, level = 0.95)[2,2] >= 1)){
        # print(paste0(name,' - ',c,'X - Slope != 1'))
      }
      confint(mod, level = 0.95)
      # print(summary(mod))
    }
    dev.off()
    colnames(OUT) <- c('covg','adj.r2','p.val','int.lo.lim','int','int.up.lim','slope.lo.lim','slope','slope.up.lim')
    write.csv(OUT, paste0('../data/bcftools_output/',d,'_true_v_callfROH_stats',name,'.csv'), row.names = FALSE)
  }
}

##### New figures accommodating multiple demo scenarios #####
##### Fig. 3 - False pos and neg at 5X and 50X #####

## False positives
text.size <- 1.25
shrink <- 2000 ## higher # here ==> narrower x-direction spread for points
alph <- 0.4
poly.alph <- 0.6
false.pos.y.lim <- 0.2
percentile <- 50
low.q <- (100-percentile)/100/2
hi.q <- 1 - (100-percentile)/100/2
diff <- 0.28 ## distance between GT and PL/PLINK
poly.wid <- 2.2 ## higher number --> narrower polygons
OUT <- NULL ## for storing means

k <- 1
while(k == 1){
pdf(paste0('../figures/all_scenarios_3_methods_false_pos_rates_5X_and_50X.pdf'), width = 10, height = 4.75)
  par(mar = c(5.1,6.1,4.1,2.1))
  plot(pl.true.v.called$covg, pl.true.v.called$false.pos.rate, ylim = c(0, false.pos.y.lim), xlab = 'Coverage', ylab = 'False positive rate', col = 'transparent', main = 'False positive rates', xlim = c(0.8,8.2), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
  axis(2, at = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), labels = c('0.0','','0.1','','0.2','','0.3'), cex.axis = text.size)
  axis(1, at = c(1:8), labels = c('5X','50X','5X','50X','5X','50X','5X','50X'), cex.axis = text.size)
  x <- 1
  for(d in c('large-1000','small','bottle','decline')){
    dat <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_vtrained_PL.csv'))
    for(c in c(5, 50)){
      temp <- dat[dat$covg == c,]
      f <- sample(c(-100:100), length(unique(temp$id)))
      f <- f/shrink+x
  
      polygon(x = c(x-diff/poly.wid, x+diff/poly.wid, x+diff/poly.wid, x-diff/poly.wid), y = c(quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(pl.col, poly.alph), border = NA)
  
      lines(x = c(x-diff/poly.wid, x+diff/poly.wid), y = c(median(temp$false.pos.rate, na.rm = TRUE), median(temp$false.pos.rate, na.rm = TRUE))
            , pch = 16, col = pl.col, lwd = 2)
      points(f, temp$false.pos.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
      OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 1))
      x <- x+1
    }
  }

  x <- 1-diff
  for(d in c('large-1000','small','bottle','decline')){
    dat <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_vtrained_GT.csv'))
    for(c in c(5, 50)){
      temp <- dat[dat$covg == c,]
      f <- sample(c(-100:100), length(unique(temp$id)))
      f <- f/shrink+x
      
      polygon(x = c(x-diff/poly.wid, x+diff/poly.wid, x+diff/poly.wid, x-diff/poly.wid), y = c(quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(gt.col, poly.alph), border = NA)
      
      lines(x = c(x-diff/poly.wid, x+diff/poly.wid), y = c(median(temp$false.pos.rate, na.rm = TRUE), median(temp$false.pos.rate, na.rm = TRUE))
            , pch = 16, col = gt.col, lwd = 2)
      points(f, temp$false.pos.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
      OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 1))
      x <- x+1
    }
  }

  x <- 1+diff
  for(d in c('large-1000','small','bottle','decline')){
    dat <- read.csv(paste0('plink_output/false_pos_neg_rates_',d,'_vtrained_plink.csv'))
    for(c in c(5, 50)){
      temp <- dat[dat$covg == c,]
      f <- sample(c(-100:100), length(unique(temp$id)))
      f <- f/shrink+x
      
      polygon(x = c(x-diff/poly.wid, x+diff/poly.wid, x+diff/poly.wid, x-diff/poly.wid), y = c(quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp$false.pos.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(plink.col, poly.alph), border = NA)
      
      lines(x = c(x-diff/poly.wid, x+diff/poly.wid), y = c(median(temp$false.pos.rate, na.rm = TRUE), median(temp$false.pos.rate, na.rm = TRUE))
            , pch = 16, col = plink.col, lwd = 2)
      points(f, temp$false.pos.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
      OUT <- rbind(OUT, c(x, mean(temp$false.pos.rate), 1))
      x <- x+1
    }
  }

  ## turn on or off vertical dashed lines between coverage levels
  abline(v = 2.5, lty = 2, col = 'grey')
  abline(v = 4.5, lty = 2, col = 'grey')
  abline(v = 6.5, lty = 2, col = 'grey')

  legend('topright', fill = c(gt.col, pl.col, plink.col), legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'),
         inset = 0.02, bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)

dev.off()
k <- k+1
}

k <- 1
while(k == 1){
  pdf(paste0('../figures/all_scenarios_3_methods_false_neg_rates_5X_and_50X.pdf'), width = 10, height = 4.75)
  par(mar = c(5.1,6.1,4.1,2.1))
  plot(pl.true.v.called$covg, pl.true.v.called$false.neg.rate, ylim = c(0, 1), xlab = 'Coverage', ylab = 'False negative rate', col = 'transparent', main = 'False negative rates', xlim = c(0.8,8.2), cex.axis = text.size, cex.lab = text.size, xaxt = 'n', yaxt = 'n')
  axis(2, at = seq(0, 1, 0.2), cex.axis = text.size)
  axis(1, at = c(1:8), labels = c('5X','50X','5X','50X','5X','50X','5X','50X'), cex.axis = text.size)
  x <- 1
  for(d in c('large-1000','small','bottle','decline')){
    dat <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_vtrained_PL.csv'))
    for(c in c(5, 50)){
      temp <- dat[dat$covg == c,]
      f <- sample(c(-100:100), length(unique(temp$id)))
      f <- f/shrink+x
      
      polygon(x = c(x-diff/poly.wid, x+diff/poly.wid, x+diff/poly.wid, x-diff/poly.wid), y = c(quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(pl.col, poly.alph), border = NA)
      
      lines(x = c(x-diff/poly.wid, x+diff/poly.wid), y = c(median(temp$false.neg.rate, na.rm = TRUE), median(temp$false.neg.rate, na.rm = TRUE))
            , pch = 16, col = pl.col, lwd = 2)
      points(f, temp$false.neg.rate, pch = 16, col = alpha(pl.col, alph), cex = 0.5) ## add individual points
      OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 1))
      x <- x+1
    }
  }
  
  x <- 1-diff
  for(d in c('large-1000','small','bottle','decline')){
    dat <- read.csv(paste0('bcftools_output/false_pos_neg_rates_',d,'_vtrained_GT.csv'))
    for(c in c(5, 50)){
      temp <- dat[dat$covg == c,]
      f <- sample(c(-100:100), length(unique(temp$id)))
      f <- f/shrink+x
      
      polygon(x = c(x-diff/poly.wid, x+diff/poly.wid, x+diff/poly.wid, x-diff/poly.wid), y = c(quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(gt.col, poly.alph), border = NA)
      
      lines(x = c(x-diff/poly.wid, x+diff/poly.wid), y = c(median(temp$false.neg.rate, na.rm = TRUE), median(temp$false.neg.rate, na.rm = TRUE))
            , pch = 16, col = gt.col, lwd = 2)
      points(f, temp$false.neg.rate, pch = 16, col = alpha(gt.col, alph), cex = 0.5) ## add individual points
      OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 1))
      x <- x+1
    }
  }
  
  x <- 1+diff
  for(d in c('large-1000','small','bottle','decline')){
    dat <- read.csv(paste0('plink_output/false_pos_neg_rates_',d,'_vtrained_plink.csv'))
    for(c in c(5, 50)){
      temp <- dat[dat$covg == c,]
      f <- sample(c(-100:100), length(unique(temp$id)))
      f <- f/shrink+x
      
      polygon(x = c(x-diff/poly.wid, x+diff/poly.wid, x+diff/poly.wid, x-diff/poly.wid), y = c(quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[1],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2],
                                                                   quantile(temp$false.neg.rate, probs = c(low.q, hi.q))[2]),
              col = alpha(plink.col, poly.alph), border = NA)
      
      lines(x = c(x-diff/poly.wid, x+diff/poly.wid), y = c(median(temp$false.neg.rate, na.rm = TRUE), median(temp$false.neg.rate, na.rm = TRUE))
            , pch = 16, col = plink.col, lwd = 2)
      points(f, temp$false.neg.rate, pch = 16, col = alpha(plink.col, alph), cex = 0.5) ## add individual points
      OUT <- rbind(OUT, c(x, mean(temp$false.neg.rate), 1))
      x <- x+1
    }
  }
  
  ## turn on or off vertical dashed lines between coverage levels
  abline(v = 2.5, lty = 2, col = 'grey')
  abline(v = 4.5, lty = 2, col = 'grey')
  abline(v = 6.5, lty = 2, col = 'grey')
  
  # legend('topright', fill = c(gt.col, pl.col, plink.col), legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK'),
  #        inset = 0.02, bg = 'gray95', box.col = 'white', border = 'transparent', cex = text.size)
  
  dev.off()
  k <- k+1
}

##### Fig. 4 - Diff between true and called f(ROH) by length bin at 15X coverage #####
pdf('/Users/Avril/Desktop/fROH_by_length_bins.pdf', width = 5, height = 4)
colors <- cbind(c('GT','PL','PLINK'), c(gt.col, pl.col, plink.col))

  for(d in c('large-1000','small','bottle','decline')){
    plot(0, 0, xlim = c(0.7, 15.3), ylim = c(-0.5, 0.25),
         col = 'transparent', main = d, xlab = '', xaxt = 'n', yaxt = 'n',
         ylab = substitute(paste('Called ',italic('F')[ROH],' - True ',italic('F')[ROH])))
    abline(h = 0, lty = 2, col = 'grey30')
    abline(v = 4, lty = 3, col = 'grey70')
    abline(v = 8, lty = 3, col = 'grey70')
    abline(v = 12, lty = 3, col = 'grey70')
    axis(1, at = c(2,6,10,14), labels = c('','','',''))
    axis(2, at = c(-0.50, -0.25, 0, 0.25))
    
    true.rohs <- read.table(paste0('slim_true_data/',d,'_true_roh_coords.txt'))
    colnames(true.rohs) <- c('id','start','end','length','true.roh.id')
    true.rohs <- true.rohs[true.rohs$length >= 100000,]
    chrom.len <- 30e6
    if(d == 'decline'){
      true.rohs <- true.rohs[true.rohs$id %notin% c(33,48),]
    }
    
    ## calculate true f(ROH) for bins
    OUT <- NULL
    for(i in unique(true.rohs$id)){
      temp <- true.rohs[true.rohs$id == i,]
      true.total <- sum(temp$length)/chrom.len
      true.b1 <- sum(temp[temp$length < b1, 'length'])/chrom.len
      true.b2 <- sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len
      true.b3 <- sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len
      true.b4 <- sum(temp[temp$length >= b3, 'length'])/chrom.len
      save <- c(i, true.total, true.b1, true.b2, true.b3, true.b4)
      OUT <- rbind(OUT, save)
    }
    true.roh.bins <- as.data.frame(OUT)
    colnames(true.roh.bins) <- c('id','true.total','true.b1','true.b2','true.b3','true.b4')
    for(c in 1:ncol(true.roh.bins)){
      true.roh.bins[,c] <- as.numeric(true.roh.bins[,c])
    }
    
    x <- 1
    
    ## bin 1
    for(m in c('GT','PL','PLINK')){
      color <- colors[colors[,1] == m, 2]
      dat <- read.csv(paste0('3_methods_results/',d,'_',m,'_frohs_by_bins.csv'))
      dat <- dat[dat$covg == 15,]
      
      DIFF <- NULL
      for(i in unique(true.roh.bins$id)){
        DIFF <- c(DIFF, dat[dat$id == i, 'bin1.froh'] - true.roh.bins[true.roh.bins$id == i, 'true.b1'])
      }
      points(x, mean(DIFF), pch = 16, col = color)
      arrows(x0 = x, x1 = x, y0 = (mean(DIFF, na.rm = TRUE) - sd(DIFF, na.rm = TRUE)),
             y1 = (mean(DIFF, na.rm = TRUE) + sd(DIFF, na.rm = TRUE)),
             lwd = 2, col = color, code=3, angle=90, length=0.08)
      if(m == 'PLINK'){
        x <- x+2
      } else{
        x <- x+1
      }
    }
    
    ## bin 2
    for(m in c('GT','PL','PLINK')){
      color <- colors[colors[,1] == m, 2]
      dat <- read.csv(paste0('3_methods_results/',d,'_',m,'_frohs_by_bins.csv'))
      dat <- dat[dat$covg == 15,]
      
      DIFF <- NULL
      for(i in unique(true.roh.bins$id)){
        DIFF <- c(DIFF, dat[dat$id == i, 'bin2.froh'] - true.roh.bins[true.roh.bins$id == i, 'true.b2'])
      }
      points(x, mean(DIFF), pch = 16, col = color)
      arrows(x0 = x, x1 = x, y0 = (mean(DIFF, na.rm = TRUE) - sd(DIFF, na.rm = TRUE)),
             y1 = (mean(DIFF, na.rm = TRUE) + sd(DIFF, na.rm = TRUE)),
             lwd = 2, col = color, code=3, angle=90, length=0.08)
      if(m == 'PLINK'){
        x <- x+2
      } else{
        x <- x+1
      }
    }
    
    ## bin 3  
    for(m in c('GT','PL','PLINK')){
      color <- colors[colors[,1] == m, 2]
      dat <- read.csv(paste0('3_methods_results/',d,'_',m,'_frohs_by_bins.csv'))
      dat <- dat[dat$covg == 15,]

      DIFF <- NULL
      for(i in unique(true.roh.bins$id)){
        DIFF <- c(DIFF, dat[dat$id == i, 'bin3.froh'] - true.roh.bins[true.roh.bins$id == i, 'true.b3'])
      }
      points(x, mean(DIFF), pch = 16, col = color)
      arrows(x0 = x, x1 = x, y0 = (mean(DIFF, na.rm = TRUE) - sd(DIFF, na.rm = TRUE)),
             y1 = (mean(DIFF, na.rm = TRUE) + sd(DIFF, na.rm = TRUE)),
             lwd = 2, col = color, code=3, angle=90, length=0.08)
      if(m == 'PLINK'){
        x <- x+2
      } else{
        x <- x+1
      }
    }
    
    ## bin 4
    for(m in c('GT','PL','PLINK')){
      color <- colors[colors[,1] == m, 2]
      dat <- read.csv(paste0('3_methods_results/',d,'_',m,'_frohs_by_bins.csv'))
      dat <- dat[dat$covg == 15,]
      
      DIFF <- NULL
      for(i in unique(true.roh.bins$id)){
        DIFF <- c(DIFF, dat[dat$id == i, 'bin4.froh'] - true.roh.bins[true.roh.bins$id == i, 'true.b4'])
      }
      points(x, mean(DIFF), pch = 16, col = color)
      arrows(x0 = x, x1 = x, y0 = (mean(DIFF, na.rm = TRUE) - sd(DIFF, na.rm = TRUE)),
             y1 = (mean(DIFF, na.rm = TRUE) + sd(DIFF, na.rm = TRUE)),
             lwd = 2, col = color, code=3, angle=90, length=0.08)
      if(m == 'PLINK'){
        x <- x+2
      } else{
        x <- x+1
      }
    }
  }
dev.off()


##### Fig. 5 - Lumping issue figure for main text #####
lims <- cbind(c('large-1000','small','bottle','decline'),
              c(9, 8, 8, 8),
              c(4.5, 5, 6, 8))

pdf('/Users/Avril/Desktop/test.pdf', width = 4.5, height = 4.5)
for(d in c('large-1000','small','bottle','decline')){
  dat <- read.csv(paste0('3_methods_results/',d,'_lumping_data_by_bins.csv'))
  
  dat[dat$method == 1, 'method'] <- 'GT'
  dat[dat$method == 2, 'method'] <- 'PL'
  dat[dat$method == 3, 'method'] <- 'PLINK'
  dat <- dat[dat$call.len >= 100e3,]
  bin.pchs <- cbind(c(1, 2, 3, 4), c(21, 22, 23, 24))
  
  pt.cex <- 1
  y.max <- as.numeric(lims[lims[,1] == d, 2])
  x.max <- as.numeric(lims[lims[,1] == d, 3])
  
  for(c in unique(dat$covg)){
    plot(0, 0, col = 'transparent', main = paste0(d,' - ',c,'X'),
         xlim = c(log(100e3), log(x.max*1e6)), xaxt = 'n',
         ylim = c(1, y.max),
         ylab = 'Number of true ROHs combined into single called ROH',
         xlab = 'Called ROH length (Mb)')
    abline(h = 1, lty = 2, col = 'grey')
    axis(1, at = log(c(100e3, 500e3, seq(0, max(dat$call.len), 1e6))), labels = c(0.1, 0.5, seq(0, max(dat$call.len), 1e6)/1e6))
    for(m in unique(dat$method)){
      sub <- dat[dat$covg == c & dat$method == m,]
      for(b in unique(sub$bin)){
        temp <- sub[sub$bin == b,]
        if(nrow(temp) > 0){
          shape <- bin.pchs[bin.pchs[,1] == b, 2]
          suppressWarnings(arrows(x0 = mean(log(temp$call.len), na.rm = TRUE), x1 = mean(log(temp$call.len), na.rm = TRUE), 
                 y0 = (mean(temp$num.true, na.rm = TRUE) - sd(temp$num.true, na.rm = TRUE)),
                 y1 = (mean(temp$num.true, na.rm = TRUE) + sd(temp$num.true, na.rm = TRUE)),
                 lwd = 1, col = temp$color[1], code=3, angle=90, length=0))
          suppressWarnings(arrows(x0 = (mean(log(temp$call.len), na.rm = TRUE) - sd(log(temp$call.len), na.rm = TRUE)), 
                 x1 = (mean(log(temp$call.len), na.rm = TRUE) + sd(log(temp$call.len), na.rm = TRUE)), 
                 y0 = mean(temp$num.true, na.rm = TRUE),
                 y1 = mean(temp$num.true, na.rm = TRUE),
                 lwd = 1, col = temp$color[1], code=3, angle=90, length=0))
          points(mean(log(temp$call.len), na.rm = TRUE), mean(temp$num.true, na.rm = TRUE), 
                 pch = shape, bg = temp$color[1], col = 'black', cex = pt.cex)
        }
      }
    }
    if(c == 5){
      legend('topleft', inset = 0.02, col = c(NA, NA, NA, 'black','black','black','black'),
             pch = c(NA, NA, NA, bin.pchs[,2]), fill = c(gt.col, pl.col, plink.col, NA, NA, NA, NA),
             border = c(NA, NA, NA, NA, NA, NA, NA), bty = 'n', pt.bg = c(NA, NA, NA, 'grey','grey','grey','grey'),
             legend = c('BCFtools Genotypes','BCFtools Likelihoods','PLINK','Short','Intermediate','Long','Very long'))
    }
  }
}
dev.off()

##### SI Fig - f(ROH) correlation calculations #####
lrg.dat <- read.csv('3_methods_results/large-1000_true_vs_called_froh_data.csv')
sml.dat <- read.csv('3_methods_results/small_true_vs_called_froh_data.csv')
bot.dat <- read.csv('3_methods_results/bottle_true_vs_called_froh_data.csv')
dec.dat <- read.csv('3_methods_results/decline_true_vs_called_froh_data.csv')

### Within each scenario, within each method, calculate among true + coverages
alph <- 0.7
pt.size <- 0.6
colors <- c(pl.col,gt.col,plink.col)

## Functions for pairs plotting
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.5)
}
onetoone_line <- function(x,y,...){
  abline(a = 0, b = 1, lty = 2, col = 'grey30')
  points(x, y, ...)
}

## Large
pdf('/Users/Avril/Desktop/large_all_covgs_corr.pdf', width = 6, height = 5.5)
col <- 1
for(m in c(4, 5, 6)){
  OUT <- NULL
  method <- gsub('.froh', '', colnames(lrg.dat)[m])
  sub <- lrg.dat[,c(1:3,m)]
  for(i in unique(sub$id)){
    temp <- sub[sub$id == i,]
    save <- c(i, temp$true.froh[1], temp[temp$covg == 5, 4],
              temp[temp$covg == 10, 4],
              temp[temp$covg == 15, 4],
              temp[temp$covg == 30, 4],
              temp[temp$covg == 50, 4])
    OUT <- rbind(OUT, save)
  }
  out <- as.data.frame(OUT)
  colnames(out) <- c('id','True','5X','10X','15X','30X','50X')
  
  pairs(out[,which(colnames(out) != 'id')], pch = 16, col = alpha(colors[col], alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size,
        gap=0, row1attop = FALSE, main = paste0('large - ',method), font.labels = 2)
  
  col <- col+1
}
dev.off()

## Small
pdf('/Users/Avril/Desktop/small_all_covgs_corr.pdf', width = 6, height = 5.5)
col <- 1
for(m in c(4, 5, 6)){
  OUT <- NULL
  method <- gsub('.froh', '', colnames(sml.dat)[m])
  sub <- sml.dat[,c(1:3,m)]
  for(i in unique(sub$id)){
    temp <- sub[sub$id == i,]
    save <- c(i, temp$true.froh[1], temp[temp$covg == 5, 4],
              temp[temp$covg == 10, 4],
              temp[temp$covg == 15, 4],
              temp[temp$covg == 30, 4],
              temp[temp$covg == 50, 4])
    OUT <- rbind(OUT, save)
  }
  out <- as.data.frame(OUT)
  colnames(out) <- c('id','True','5X','10X','15X','30X','50X')
  
  pairs(out[,which(colnames(out) != 'id')], pch = 16, col = alpha(colors[col], alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size,
        gap=0, row1attop = FALSE, main = paste0('small - ',method), font.labels = 2)
  
  col <- col+1
}
dev.off()

## Bottle
pdf('/Users/Avril/Desktop/bottle_all_covgs_corr.pdf', width = 6, height = 5.5)
col <- 1
for(m in c(4, 5, 6)){
  OUT <- NULL
  method <- gsub('.froh', '', colnames(bot.dat)[m])
  sub <- bot.dat[,c(1:3,m)]
  for(i in unique(sub$id)){
    temp <- sub[sub$id == i,]
    save <- c(i, temp$true.froh[1], temp[temp$covg == 5, 4],
              temp[temp$covg == 10, 4],
              temp[temp$covg == 15, 4],
              temp[temp$covg == 30, 4],
              temp[temp$covg == 50, 4])
    OUT <- rbind(OUT, save)
  }
  out <- as.data.frame(OUT)
  colnames(out) <- c('id','True','5X','10X','15X','30X','50X')
  
  pairs(out[,which(colnames(out) != 'id')], pch = 16, col = alpha(colors[col], alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size,
        gap=0, row1attop = FALSE, main = paste0('bottle - ',method), font.labels = 2)
  
  col <- col+1
}
dev.off()

## Decline
pdf('/Users/Avril/Desktop/decline_all_covgs_corr.pdf', width = 6, height = 5.5)
col <- 1
for(m in c(4, 5, 6)){
  OUT <- NULL
  method <- gsub('.froh', '', colnames(dec.dat)[m])
  sub <- dec.dat[,c(1:3,m)]
  for(i in unique(sub$id)){
    temp <- sub[sub$id == i,]
    save <- c(i, temp$true.froh[1], temp[temp$covg == 5, 4],
              temp[temp$covg == 10, 4],
              temp[temp$covg == 15, 4],
              temp[temp$covg == 30, 4],
              temp[temp$covg == 50, 4])
    OUT <- rbind(OUT, save)
  }
  out <- as.data.frame(OUT)
  colnames(out) <- c('id','True','5X','10X','15X','30X','50X')
  
  pairs(out[,which(colnames(out) != 'id')], pch = 16, col = alpha(colors[col], alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size,
        gap=0, row1attop = FALSE, main = paste0('decline - ',method), font.labels = 2)
  
  col <- col+1
}
dev.off()

## Large
pdf('/Users/Avril/Desktop/large_all_methods_corr.pdf', width = 6, height = 5.5)
col <- 1
for(c in c(5, 10, 15, 30, 50)){
  sub <- lrg.dat[lrg.dat$covg == c,]
  
  pairs(sub[,which(colnames(sub) %notin% c('id','covg'))], pch = 16, col = alpha('springgreen4', alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size, cex.labels = 1.75,
        gap=0, row1attop = FALSE, main = paste0('large - ',c,'X'), font.labels = 2,
        labels = c('True','BCFtools\nLikelihoods','BCFtools\nGenotypes','PLINK'))
}
dev.off()

## Small
pdf('/Users/Avril/Desktop/small_all_methods_corr.pdf', width = 6, height = 5.5)
col <- 1
for(c in c(5, 10, 15, 30, 50)){
  sub <- sml.dat[sml.dat$covg == c,]
  
  pairs(sub[,which(colnames(sub) %notin% c('id','covg'))], pch = 16, col = alpha('springgreen4', alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size, cex.labels = 1.75,
        gap=0, row1attop = FALSE, main = paste0('large - ',c,'X'), font.labels = 2,
        labels = c('True','BCFtools\nLikelihoods','BCFtools\nGenotypes','PLINK'))
}
dev.off()

## Bottle
pdf('/Users/Avril/Desktop/bottle_all_methods_corr.pdf', width = 6, height = 5.5)
for(c in c(5, 10, 15, 30, 50)){
  sub <- bot.dat[bot.dat$covg == c,]
  
  pairs(sub[,which(colnames(sub) %notin% c('id','covg'))], pch = 16, col = alpha('springgreen4', alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size, cex.labels = 1.75,
        gap=0, row1attop = FALSE, main = paste0('large - ',c,'X'), font.labels = 2,
        labels = c('True','BCFtools\nLikelihoods','BCFtools\nGenotypes','PLINK'))
}
dev.off()

## Decline
pdf('/Users/Avril/Desktop/decline_all_methods_corr.pdf', width = 6, height = 5.5)
col <- 1
for(c in c(5, 10, 15, 30, 50)){
  sub <- dec.dat[dec.dat$covg == c,]
  
  pairs(sub[,which(colnames(sub) %notin% c('id','covg'))], pch = 16, col = alpha('springgreen4', alph), 
        upper.panel = panel.cor, lower.panel = onetoone_line, cex = pt.size, cex.labels = 1.75,
        gap=0, row1attop = FALSE, main = paste0('large - ',c,'X'), font.labels = 2,
        labels = c('True','BCFtools\nLikelihoods','BCFtools\nGenotypes','PLINK'))
}
dev.off()
