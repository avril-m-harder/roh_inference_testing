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
  scen.name <- demo.names[demo.names[,1] == d, 2]
  for(h in c('vtrained')){
    print(paste0(d,' - ',h))
    
    ##### 1A. Read in data and summary files #####
    ### True ROH information
    true.rohs <- read.table(paste0('slim_true_data/',d,'_true_roh_coords.txt'))
    colnames(true.rohs) <- c('id','start','end','length','true.roh.id')
    true.rohs <- true.rohs[true.rohs$length >= 100000,]
    chrom.len <- 30e6
    if(d == 'decline'){
      true.rohs <- true.rohs[true.rohs$id %notin% c(33,48),]
    }

    ##### >>> HWE results #####
    ### Read in ROH calling results
    hwe.bcf.res <- read.table('bcftools_output/hwe_filtered/bcftoolsROH_all_coordinates.txt',
                          header = TRUE, sep = '\t')
    hwe.bcf.res$covg <- as.numeric(gsub('x', '', hwe.bcf.res$covg))
    hwe.bcf.res$id <- gsub('i', '', hwe.bcf.res$id)
    hwe.bcf.res <- hwe.bcf.res[hwe.bcf.res$demo == d,]
    if(d == 'decline'){
      hwe.bcf.res <- hwe.bcf.res[hwe.bcf.res$id %notin% c(33,48),]
    }

    ## GT
    hwe.bcf.gt.res <- hwe.bcf.res[hwe.bcf.res$method == 'GT' & hwe.bcf.res$hmm == h,]
    hwe.bcf.gt.res <- hwe.bcf.gt.res[hwe.bcf.gt.res$length >= 100000,] ## applying 100kb filter
    hwe.bcf.gt.ids <- unique(hwe.bcf.gt.res$id)

    ## PL
    hwe.bcf.pl.res <- hwe.bcf.res[hwe.bcf.res$method == 'PL' & hwe.bcf.res$hmm == h,]
    hwe.bcf.pl.res <- hwe.bcf.pl.res[hwe.bcf.pl.res$length >= 100000,] ## applying 100kb filter
    hwe.bcf.pl.ids <- unique(hwe.bcf.pl.res$id)

    ## PLINK
    hwe.plink.res <- read.table(paste0('plink_output/hwe_round1/',d,'_PLINK_all_coordinates.txt'), header = TRUE)
    ##### !!! Update this for all demo scenarios, covgs if necessary !!! #####
    ## final selection = default settings
    hwe.plink.res <- hwe.plink.res[hwe.plink.res$phwh == 1 & hwe.plink.res$phwm == 5 & hwe.plink.res$phws == 50 &
                             hwe.plink.res$phwt == 0.05 & hwe.plink.res$phzs == 100 & hwe.plink.res$phzg == 1000,]
    hwe.plink.res$length <- hwe.plink.res$end - hwe.plink.res$start + 1
    hwe.plink.res$id <- gsub('i', '', hwe.plink.res$id)
    if(d == 'decline'){
      hwe.plink.res <- hwe.plink.res[hwe.plink.res$id %notin% c(33,48),]
    }
    hwe.plink.ids <- unique(hwe.plink.res$id)
    hwe.plink.res$covg <- as.numeric(gsub('x', '', hwe.plink.res$covg))

    ### Write/read in called vs. true f(ROH) results
    OUT <- NULL
    for(i in unique(true.rohs$id)){
      true.froh <- sum(true.rohs[true.rohs$id == i, 'length'])/chrom.len
      for(c in unique(hwe.bcf.gt.res$covg)){
        gt.froh <- sum(hwe.bcf.gt.res[hwe.bcf.gt.res$id == i & hwe.bcf.gt.res$covg == c, 'length'])/chrom.len
        pl.froh <- sum(hwe.bcf.pl.res[hwe.bcf.pl.res$id == i & hwe.bcf.pl.res$covg == c, 'length'])/chrom.len
        pk.froh <- sum(hwe.plink.res[hwe.plink.res$id == i & hwe.plink.res$covg == c, 'length'])/chrom.len
        save <- c(i, c, true.froh, pl.froh, gt.froh, pk.froh)
        OUT <- rbind(OUT, save)
      }
    }
    colnames(OUT) <- c('id','covg','true.froh','pl.froh','gt.froh','plink.froh')
    write.csv(OUT, paste0('3_methods_results/hwe_filtering_',d,'_true_vs_called_froh_data.csv'), row.names = FALSE)
    hwe.froh.stats <- read.csv(paste0('3_methods_results/hwe_filtering_',d,'_true_vs_called_froh_data.csv'), header = TRUE)
    
    ##### >>> non-HWE results #####
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
    bcf.gt.ids <- unique(bcf.gt.res$id)
    
    
    ## PL
    bcf.pl.res <- bcf.res[bcf.res$method == 'PL' & bcf.res$hmm == h,]
    bcf.pl.res <- bcf.pl.res[bcf.pl.res$length >= 100000,] ## applying 100kb filter
    bcf.pl.ids <- unique(bcf.pl.res$id)
    
    ## PLINK
    plink.res <- read.table(paste0('plink_output/round1/',d,'_PLINK_all_coordinates.txt'), header = TRUE)
    ##### !!! Update this for all demo scenarios, covgs if necessary !!! #####
    ## final selection = default settings
    plink.res <- plink.res[plink.res$phwh == 1 & plink.res$phwm == 5 & plink.res$phws == 50 &
                                     plink.res$phwt == 0.05 & plink.res$phzs == 100 & plink.res$phzg == 1000,]
    plink.res$length <- plink.res$end - plink.res$start + 1
    plink.res$id <- gsub('i', '', plink.res$id)
    if(d == 'decline'){
      plink.res <- plink.res[plink.res$id %notin% c(33,48),]
    }
    plink.ids <- unique(plink.res$id)
    plink.res$covg <- as.numeric(gsub('x', '', plink.res$covg))
    
    ## Get list of IDs with results across all methods/approaches and apply to all sets
    ids <- table(c(hwe.bcf.gt.ids, hwe.bcf.pl.ids, hwe.plink.ids,
                 bcf.gt.ids, bcf.pl.ids, plink.ids))
    ids <- names(ids[which(ids == 6)])
    
    bcf.gt.res <- bcf.gt.res[which(bcf.gt.res$id %in% ids),]
    bcf.pl.res <- bcf.pl.res[which(bcf.pl.res$id %in% ids),]
    plink.res <- plink.res[which(plink.res$id %in% ids),]
    hwe.bcf.gt.res <- hwe.bcf.gt.res[which(hwe.bcf.gt.res$id %in% ids),]
    hwe.bcf.pl.res <- hwe.bcf.pl.res[which(hwe.bcf.pl.res$id %in% ids),]
    hwe.plink.res <- hwe.plink.res[which(hwe.plink.res$id %in% ids),]
    gt.frohs <- gt.frohs[which(gt.frohs$id %in% ids),]
    pl.frohs <- pl.frohs[which(pl.frohs$id %in% ids),]
    plink.frohs <- plink.frohs[which(plink.frohs$id %in% ids),]
    hwe.gt.frohs <- hwe.gt.frohs[which(hwe.gt.frohs$id %in% ids),]
    hwe.pl.frohs <- hwe.pl.frohs[which(hwe.pl.frohs$id %in% ids),]
    hwe.plink.frohs <- hwe.plink.frohs[which(hwe.plink.frohs$id %in% ids),]
    true.frohs <- true.frohs[which(true.frohs$id %in% ids),]
    
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

    ##### 94A. Plotting f(ROH) by length bins #####
    b1 <- 500e3
    b2 <- 1e6
    b3 <- 2e6

    ymin <- 0
    ymax <- 1
    alph <- 0.8

    k <- 1
    while(k == 1){
      pdf(paste0('../figures/hwe_filtering/',d,'/',d,'_',h,'_fROH_by_length_bins_HWE_vs_nonHWE.pdf'), width = 12, height = 3.25)
      par(mfrow = c(1,4))

      ## GT
      ## HWE
      OUT <- NULL
      for(i in unique(hwe.bcf.gt.res$id)){
        for(c in unique(hwe.bcf.gt.res$covg)){
          temp <- hwe.bcf.gt.res[hwe.bcf.gt.res$id == i & hwe.bcf.gt.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      hwe.gt.frohs <- as.data.frame(OUT)
      colnames(hwe.gt.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(hwe.gt.frohs)){
        hwe.gt.frohs[,c] <- as.numeric(hwe.gt.frohs[,c])
      }
      
      ## non-HWE
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

      ## GT plots
      hwe.gt.frohs <- hwe.gt.frohs[order(hwe.gt.frohs$id),]
      gt.frohs <- gt.frohs[order(gt.frohs$id),]
      if(all(hwe.gt.frohs$id != gt.frohs$id)){
        print(paste0('GT ID issues'))
      }
      
      for(c in unique(gt.frohs$covg)){
        for(b in c('bin1','bin2','bin3','bin4')){
          col <- grep(paste0(b,'.froh'), colnames(gt.frohs))
          hwe.col <- grep(paste0(b,'.froh'), colnames(hwe.gt.frohs))
          
          plot(gt.frohs[,col], hwe.gt.frohs[,hwe.col], main = paste0(b,'\n',c,'X'), 
               xlab = substitute(paste('Non-HWE called ',italic('F')[ROH])),
               ylab = substitute(paste('HWE called ',italic('F')[ROH])), col = 'transparent')
            abline(0, 1, lty = 2, col = 'darkgrey')
            points(gt.frohs[gt.frohs$covg == c, col], hwe.gt.frohs[hwe.gt.frohs$covg == c, hwe.col],
                   pch = 16, col = alpha(gt.col, alph))
        }
      }


      ## PL
      ## HWE
      OUT <- NULL
      for(i in unique(hwe.bcf.pl.res$id)){
        for(c in unique(hwe.bcf.pl.res$covg)){
          temp <- hwe.bcf.pl.res[hwe.bcf.pl.res$id == i & hwe.bcf.pl.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      hwe.pl.frohs <- as.data.frame(OUT)
      colnames(hwe.pl.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(hwe.pl.frohs)){
        hwe.pl.frohs[,c] <- as.numeric(hwe.pl.frohs[,c])
      }
      
      ## non-HWE
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
      
      ## PL plots
      hwe.pl.frohs <- hwe.pl.frohs[order(hwe.pl.frohs$id),]
      pl.frohs <- pl.frohs[order(pl.frohs$id),]
      if(all(hwe.pl.frohs$id != pl.frohs$id)){
        print(paste0('PL ID issues'))
      }
      
      for(c in unique(pl.frohs$covg)){
        for(b in c('bin1','bin2','bin3','bin4')){
          col <- grep(paste0(b,'.froh'), colnames(pl.frohs))
          hwe.col <- grep(paste0(b,'.froh'), colnames(hwe.pl.frohs))
          
          plot(pl.frohs[,col], hwe.pl.frohs[,hwe.col], main = paste0(b,'\n',c,'X'), 
               xlab = substitute(paste('Non-HWE called ',italic('F')[ROH])),
               ylab = substitute(paste('HWE called ',italic('F')[ROH])), col = 'transparent')
          abline(0, 1, lty = 2, col = 'darkgrey')
          points(pl.frohs[pl.frohs$covg == c, col], hwe.pl.frohs[hwe.pl.frohs$covg == c, hwe.col],
                 pch = 16, col = alpha(pl.col, alph))
        }
      }

      ## PLINK
      ## HWE
      OUT <- NULL
      for(i in unique(hwe.plink.res$id)){
        for(c in unique(hwe.plink.res$covg)){
          temp <- hwe.plink.res[hwe.plink.res$id == i & hwe.plink.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      hwe.plink.frohs <- as.data.frame(OUT)
      colnames(hwe.plink.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(hwe.plink.frohs)){
        hwe.plink.frohs[,c] <- as.numeric(hwe.plink.frohs[,c])
      }
      
      ## non-HWE
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
      
      ## PLINK plots
      hwe.plink.frohs <- hwe.plink.frohs[order(hwe.plink.frohs$id),]
      plink.frohs <- plink.frohs[order(plink.frohs$id),]
      if(all(hwe.plink.frohs$id != plink.frohs$id)){
        print(paste0('PLINK ID issues'))
      }
      
      for(c in unique(plink.frohs$covg)){
        for(b in c('bin1','bin2','bin3','bin4')){
          col <- grep(paste0(b,'.froh'), colnames(plink.frohs))
          hwe.col <- grep(paste0(b,'.froh'), colnames(hwe.plink.frohs))
          
          plot(plink.frohs[,col], hwe.plink.frohs[,hwe.col], main = paste0(b,'\n',c,'X'), 
               xlab = substitute(paste('Non-HWE called ',italic('F')[ROH])),
               ylab = substitute(paste('HWE called ',italic('F')[ROH])), col = 'transparent')
          abline(0, 1, lty = 2, col = 'darkgrey')
          points(plink.frohs[plink.frohs$covg == c, col], hwe.plink.frohs[hwe.plink.frohs$covg == c, hwe.col],
                 pch = 16, col = alpha(plink.col, alph))
        }
      }

      k <- k+1
      dev.off()
    }
    
    k <- 1
    while(k == 1){
      pdf(paste0('../figures/hwe_filtering/',d,'/',d,'_',h,'_fROH_by_length_bins_HWE_vs_nonHWE_wTruefROH.pdf'), width = 12, height = 3.25)
      par(mfrow = c(1,4))
      
      ## True
      OUT <- NULL
      for(i in unique(true.rohs$id)){
        temp <- true.rohs[true.rohs$id == i,]
        save <- c(i, sum(temp$length)/chrom.len, nrow(temp),
                  sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                  sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                  sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                  sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
        OUT <- rbind(OUT, save)
      }
      true.frohs <- as.data.frame(OUT)
      colnames(true.frohs) <- c('id','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(true.frohs)){
        true.frohs[,c] <- as.numeric(true.frohs[,c])
      }
      true.frohs <- true.frohs[true.frohs$id %in% ids,]
      
      ## GT
      ## HWE
      OUT <- NULL
      for(i in unique(hwe.bcf.gt.res$id)){
        for(c in unique(hwe.bcf.gt.res$covg)){
          temp <- hwe.bcf.gt.res[hwe.bcf.gt.res$id == i & hwe.bcf.gt.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      hwe.gt.frohs <- as.data.frame(OUT)
      colnames(hwe.gt.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(hwe.gt.frohs)){
        hwe.gt.frohs[,c] <- as.numeric(hwe.gt.frohs[,c])
      }
      
      ## non-HWE
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
      
      ## GT plots
      hwe.gt.frohs <- hwe.gt.frohs[order(hwe.gt.frohs$id),]
      gt.frohs <- gt.frohs[order(gt.frohs$id),]
      if(all(hwe.gt.frohs$id != gt.frohs$id) | all(unique(hwe.gt.frohs$id) != true.frohs$id)){
        print(paste0('GT ID issues'))
      }
      
      for(c in unique(gt.frohs$covg)){
        for(b in c('bin1','bin2','bin3','bin4')){
          col <- grep(paste0(b,'.froh'), colnames(gt.frohs))
          hwe.col <- grep(paste0(b,'.froh'), colnames(hwe.gt.frohs))
          true.col <- grep(paste0(b,'.froh'), colnames(true.frohs))
          
          plot(0, 0, xlim = c(min(true.frohs[,true.col], hwe.gt.frohs[,hwe.col], gt.frohs[,col]), 
                              max(true.frohs[,true.col], hwe.gt.frohs[,hwe.col], gt.frohs[,col])),
               ylim = c(min(true.frohs[,true.col], hwe.gt.frohs[,hwe.col], gt.frohs[,col]), 
                        max(true.frohs[,true.col], hwe.gt.frohs[,hwe.col], gt.frohs[,col])),
               main = paste0(b,'\n',c,'X'), 
               xlab = substitute(paste('True ',italic('F')[ROH])),
               ylab = substitute(paste('Called ',italic('F')[ROH])), col = 'transparent')
          abline(0, 1, lty = 2, col = 'darkgrey')
          points(true.frohs[, true.col], hwe.gt.frohs[hwe.gt.frohs$covg == c, hwe.col],
                 pch = 16, col = alpha(gt.col, alph))
          abline(lm(hwe.gt.frohs[hwe.gt.frohs$covg == c, hwe.col] ~ true.frohs[, true.col]), col = gt.col, lty = 3)
          points(true.frohs[, true.col], gt.frohs[gt.frohs$covg == c, col],
                 pch = 17, col = alpha(gt.col, alph))
          abline(lm(gt.frohs[gt.frohs$covg == c, col] ~ true.frohs[, true.col]), col = gt.col)
          if(b == 'bin1'){
            par(xpd = TRUE)
            legend('topleft', inset = c(-0.1, -0.2), col = gt.col, lty = c(1,3), legend = c('non-HWE','HWE'), bty = 'n')
            par(xpd = FALSE)
          }
        }
      }
      
      
      ## PL
      ## HWE
      OUT <- NULL
      for(i in unique(hwe.bcf.pl.res$id)){
        for(c in unique(hwe.bcf.pl.res$covg)){
          temp <- hwe.bcf.pl.res[hwe.bcf.pl.res$id == i & hwe.bcf.pl.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      hwe.pl.frohs <- as.data.frame(OUT)
      colnames(hwe.pl.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(hwe.pl.frohs)){
        hwe.pl.frohs[,c] <- as.numeric(hwe.pl.frohs[,c])
      }
      
      ## non-HWE
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
      
      ## PL plots
      hwe.pl.frohs <- hwe.pl.frohs[order(hwe.pl.frohs$id),]
      pl.frohs <- pl.frohs[order(pl.frohs$id),]
      if(all(hwe.pl.frohs$id != pl.frohs$id) | all(unique(hwe.pl.frohs$id) != true.frohs$id)){
        print(paste0('PL ID issues'))
      }
      
      for(c in unique(pl.frohs$covg)){
        for(b in c('bin1','bin2','bin3','bin4')){
          col <- grep(paste0(b,'.froh'), colnames(pl.frohs))
          hwe.col <- grep(paste0(b,'.froh'), colnames(hwe.pl.frohs))
          true.col <- grep(paste0(b,'.froh'), colnames(true.frohs))
          
          plot(0, 0, xlim = c(min(true.frohs[,true.col], hwe.pl.frohs[,hwe.col], pl.frohs[,col]), 
                              max(true.frohs[,true.col], hwe.pl.frohs[,hwe.col], pl.frohs[,col])),
               ylim = c(min(true.frohs[,true.col], hwe.pl.frohs[,hwe.col], pl.frohs[,col]), 
                        max(true.frohs[,true.col], hwe.pl.frohs[,hwe.col], pl.frohs[,col])),
               main = paste0(b,'\n',c,'X'), 
               xlab = substitute(paste('True ',italic('F')[ROH])),
               ylab = substitute(paste('Called ',italic('F')[ROH])), col = 'transparent')
          abline(0, 1, lty = 2, col = 'darkgrey')
          points(true.frohs[, true.col], hwe.pl.frohs[hwe.pl.frohs$covg == c, hwe.col],
                 pch = 16, col = alpha(pl.col, alph))
          abline(lm(hwe.pl.frohs[hwe.pl.frohs$covg == c, hwe.col] ~ true.frohs[, true.col]), col = pl.col, lty = 3)
          points(true.frohs[, true.col], pl.frohs[pl.frohs$covg == c, col],
                 pch = 17, col = alpha(pl.col, alph))
          abline(lm(pl.frohs[pl.frohs$covg == c, col] ~ true.frohs[, true.col]), col = pl.col)
          if(b == 'bin1'){
            par(xpd = TRUE)
            legend('topleft', inset = c(-0.1, -0.2), col = pl.col, lty = c(1,3), legend = c('non-HWE','HWE'), bty = 'n')
            par(xpd = FALSE)
          }
        }
      }
      
      ## PLINK
      ## HWE
      OUT <- NULL
      for(i in unique(hwe.plink.res$id)){
        for(c in unique(hwe.plink.res$covg)){
          temp <- hwe.plink.res[hwe.plink.res$id == i & hwe.plink.res$covg == c,]
          save <- c(i, c, sum(temp$length)/chrom.len, nrow(temp),
                    sum(temp[temp$length < b1, 'length'])/chrom.len, nrow(temp[temp$length < b1,]),
                    sum(temp[temp$length >= b1 & temp$length < b2, 'length'])/chrom.len, nrow(temp[temp$length >= b1 & temp$length < b2,]),
                    sum(temp[temp$length >= b2 & temp$length < b3, 'length'])/chrom.len, nrow(temp[temp$length >= b2 & temp$length < b3,]),
                    sum(temp[temp$length >= b3, 'length'])/chrom.len, nrow(temp[temp$length >= b3,]))
          OUT <- rbind(OUT, save)
        }
      }
      hwe.plink.frohs <- as.data.frame(OUT)
      colnames(hwe.plink.frohs) <- c('id','covg','froh','n.rohs','bin1.froh','bin1.n','bin2.froh','bin2.n','bin3.froh','bin3.n','bin4.froh','bin4.n')
      for(c in 2:ncol(hwe.plink.frohs)){
        hwe.plink.frohs[,c] <- as.numeric(hwe.plink.frohs[,c])
      }
      
      ## non-HWE
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
      
      ## PLINK plots
      hwe.plink.frohs <- hwe.plink.frohs[order(hwe.plink.frohs$id),]
      plink.frohs <- plink.frohs[order(plink.frohs$id),]
      if(all(hwe.plink.frohs$id != plink.frohs$id) | all(unique(hwe.plink.frohs$id) != true.frohs$id)){
        print(paste0('PLINK ID issues'))
      }
      
      for(c in unique(plink.frohs$covg)){
        for(b in c('bin1','bin2','bin3','bin4')){
          col <- grep(paste0(b,'.froh'), colnames(plink.frohs))
          hwe.col <- grep(paste0(b,'.froh'), colnames(hwe.plink.frohs))
          true.col <- grep(paste0(b,'.froh'), colnames(true.frohs))
          
          plot(0, 0, xlim = c(min(true.frohs[,true.col], hwe.plink.frohs[,hwe.col], plink.frohs[,col]), 
                              max(true.frohs[,true.col], hwe.plink.frohs[,hwe.col], plink.frohs[,col])),
               ylim = c(min(true.frohs[,true.col], hwe.plink.frohs[,hwe.col], plink.frohs[,col]), 
                        max(true.frohs[,true.col], hwe.plink.frohs[,hwe.col], plink.frohs[,col])),
               main = paste0(b,'\n',c,'X'), 
               xlab = substitute(paste('True ',italic('F')[ROH])),
               ylab = substitute(paste('Called ',italic('F')[ROH])), col = 'transparent')
          abline(0, 1, lty = 2, col = 'darkgrey')
          points(true.frohs[, true.col], hwe.plink.frohs[hwe.plink.frohs$covg == c, hwe.col],
                 pch = 16, col = alpha(plink.col, alph))
          abline(lm(hwe.plink.frohs[hwe.plink.frohs$covg == c, hwe.col] ~ true.frohs[, true.col]), col = plink.col, lty = 3)
          points(true.frohs[, true.col], plink.frohs[plink.frohs$covg == c, col],
                 pch = 17, col = alpha(plink.col, alph))
          abline(lm(plink.frohs[plink.frohs$covg == c, col] ~ true.frohs[, true.col]), col = plink.col)
          if(b == 'bin1'){
            par(xpd = TRUE)
            legend('topleft', inset = c(-0.1, -0.2), col = plink.col, lty = c(1,3), legend = c('non-HWE','HWE'), bty = 'n')
            par(xpd = FALSE)
          }
        }
      }
      
      k <- k+1
      dev.off()
    }
  }
}
