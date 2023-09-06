### Script to calculate and visualize heterozygosity in windows alongside true and called ROHs

library(scales)
library(ghibli)
library(vcfR)
`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/')

### set colors
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
coverages <- cbind(c(5, 10, 15, 30, 50),
                   c('05x','10x','15x','30x','50x'))

## select Viterbi-trained data
h <- 'vtrained'

## set window size
wind.size <- 1e6

true.fns <- list.files(path = 'data/slim_true_data/', pattern = 'true_roh_coords')
demos <- do.call(rbind, strsplit(true.fns, split = '_'))[,1]

for(d in demos){
  dir.create(paste0('../figures/',d), showWarnings = FALSE)
  scen.name <- demo.names[demo.names[,1] == d, 2]
  
  ##### 1. Read in data and summary files #####
  ### True ROH information
  true.rohs <- read.table(paste0('data/slim_true_data/',d,'_true_roh_coords.txt'))
  colnames(true.rohs) <- c('id','start','end','length','true.roh.id')
  true.rohs <- true.rohs[true.rohs$length >= 100000,]
  chrom.len <- 30e6
  if(d == 'decline'){
    true.rohs <- true.rohs[true.rohs$id %notin% c(33,48),]
  }
  
  ### Read in ROH calling results
  bcf.res <- read.table('data/bcftools_output/bcftoolsROH_all_coordinates.txt',
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
  plink.res <- read.table(paste0('data/plink_final_iteration/',d,'_PLINK_all_coordinates.txt'), header = TRUE)
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
  
  
  ##### 2. Looping over coverage levels, read in VCF and calculate heterozygosity #####
  for(cov in 1:nrow(coverages)){
    
    vcf <- read.vcfR(paste0('data/filtered_snp_vcfs/',d,'_sample_cvg_',coverages[cov,2],'_filtered_SNPs.vcf'))
    gts <- extract.gt(vcf)
    colnames(gts) <- gsub('i', '', colnames(gts))
    pos <- vcf@fix[,2]
    
    ##### 99. Plotting individual true ROHs and called ROHs for all 3 analyses and coverage levels #####
    wid <- 3 ## line widths in plot ~ or ~
    hit <- 0.1 ## polygon height
    pt.cex <- 0.4
    lwd <- 0.5
    cov.name <- coverages[cov,2]
    
    pdf(paste0('/figures/',d,'/',d,'_',h,'_',cov.name,'_individual_true_and_called_ROHs_across_chromosome.pdf'), 
        width = 15, height = 9)
    par(mar = c(5.1, 10.1, 4.1, 2.1))
    for(i in sort(as.numeric(unique(bcf.pl.res$id)))){
      id.gts <- cbind(pos, gts[,which(colnames(gts) == i)])
      id.gts[id.gts[,2] == '0/0', 2] <- 0
      id.gts[id.gts[,2] == '0|0', 2] <- 0
      id.gts[id.gts[,2] == '1/1', 2] <- 0
      id.gts[id.gts[,2] == '1|1', 2] <- 0
      id.gts[id.gts[,2] == '0/1', 2] <- 1
      id.gts[id.gts[,2] == '1/0', 2] <- 1
      id.gts[id.gts[,2] == '0|1', 2] <- 1
      id.gts[id.gts[,2] == '1|0', 2] <- 1
      
      ## calculate window heterozygosity values
      s <- 1
      e <- s + wind.size - 1
      OUT <- NULL
      while(e <= chrom.len){
        mid <- e - (wind.size/2)
        het <- nrow(id.gts[which(as.numeric(id.gts[,1]) >= s &
                                 as.numeric(id.gts[,1]) <= e &
                                 as.numeric(id.gts[,2] == 1)),])
        if(is.null(het)){
          het <- 0
        }
        all <- nrow(id.gts[which(as.numeric(id.gts[,1]) >= s &
                                   as.numeric(id.gts[,1]) <= e),])
        save <- c(mid, (het/all))
        OUT <- rbind(OUT, save)
        s <- s + wind.size
        e <- s + wind.size - 1
      }
            
      sub.pl <- bcf.pl.res[bcf.pl.res$id == i & bcf.pl.res$covg == coverages[cov,1],]
      sub.gt <- bcf.gt.res[bcf.gt.res$id == i & bcf.gt.res$covg == coverages[cov,1],]
      sub.pk <- plink.res[plink.res$id == i & plink.res$covg == coverages[cov,1],]
      sub.true <- true.rohs[true.rohs$id == i,]
      
      plot(0,0, xlim = c(1, chrom.len), ylim = c(0.85, 5.65), col = 'transparent', yaxt = 'n', xlab = 'Chromosome position (bp)', 
           ylab = '', main = '', bty = 'n')
      lines(x = c(-2e6, 40e6), y = c(1.5, 1.5), lty = 2)
      # lines(x = c(-2e6, 40e6), y = c(2.5, 2.5), lty = 2)
      # lines(x = c(-2e6, 40e6), y = c(3.5, 3.5), lty = 2)
      axis(2, at = c(1:4), labels = c('True','PLINK','BCFtools\nLikelihoods','BCFtools\nGenotypes'), las = 2)
      axis(2, at = c(4.5, 5, 5.5), labels = c('0','0.5','1'), las = 2)
      mtext(text = 'Heterozygosity', side = 2, las= 2, at = 5, line = 3)
      for(r in 1:nrow(sub.true)){
        # lines(x = c(sub.true$start[r], sub.true$end[r]), y = c(1, 1), col = 'black', lwd = wid)
        polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1, 1, 18, 18),
                col = alpha('lightgrey', 0.3), border = NA)
        polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1-hit, 1-hit, 1+hit, 1+hit),
                col = 'black', border = NA)
      }
      y <- 2
      for(r in 1:nrow(sub.pk)){
        polygon(x = c(sub.pk$start[r], sub.pk$end[r], sub.pk$end[r], sub.pk$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = plink.col, border = NA)
      }
      y <- y+1
      for(r in 1:nrow(sub.pl)){
        polygon(x = c(sub.pl$start[r], sub.pl$end[r], sub.pl$end[r], sub.pl$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = pl.col, border = NA)
      }
      y <- y+1
      for(r in 1:nrow(sub.gt)){
        polygon(x = c(sub.gt$start[r], sub.gt$end[r], sub.gt$end[r], sub.gt$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
                col = gt.col, border = NA)
      }

      OUT[,2] <- OUT[,2] + 4.5
      points(OUT[,1], OUT[,2], pch = 16, cex = pt.cex, col = 'springgreen4')
      lines(OUT[,1], OUT[,2], lwd = lwd, col = 'springgreen4')
      lines(x = c(1,chrom.len), y = c(4.5, 4.5))
      
      # mtext('PLINK', side = 2, line = 3, at = 14)
      # mtext('BCFtools\nGenotypes', side = 2, line = 3, at = 9)
      # mtext('BCFtools\nLikelihoods', side = 2, line = 3, at = 4)
    }
    dev.off()
    
    ##### >>> 99A. Plotting individual true ROHs and called ROHs for all 3 analyses and coverage levels - specific individual and window #####
    #     wid <- 3 ## line widths in plot ~ or ~
    #     hit <- 0.1 ## polygon height
    #     xmin <- 1.5e7
    #     xmax <- 2.5e7
    #
    #     # pdf(paste0('../figures/',d,'/',figure.ct,'_',d,'_',h,'_individual_true_and_called_ROHs_across_chromosome.pdf'), width = 15, height = 7)
    #     par(mar = c(5.1, 5.1, 4.1, 2.1))
    #     for(i in unique(true.rohs$id)){
    #       sub.pl <- bcf.pl.res[bcf.pl.res$id == i,]
    #       sub.gt <- bcf.gt.res[bcf.gt.res$id == i,]
    #       sub.pk <- plink.res[plink.res$id == i,]
    #       sub.true <- true.rohs[true.rohs$id == i,]
    #       plot(0,0, ylim = c(1, 16), col = 'transparent', yaxt = 'n', xlab = 'Chromosome position (bp)', ylab = '', main = i,
    #            xlim = c(xmin, xmax))
    #       lines(x = c(-2e6, 40e6), y = c(1.5, 1.5), lty = 2)
    #       lines(x = c(-2e6, 40e6), y = c(6.5, 6.5), lty = 2)
    #       lines(x = c(-2e6, 40e6), y = c(11.5, 11.5), lty = 2)
    #       axis(2, at = c(1:16), labels = c('True','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X','5X','10X','15X','30X','50X'), las = 2)
    #       for(r in 1:nrow(sub.true)){
    #         # lines(x = c(sub.true$start[r], sub.true$end[r]), y = c(1, 1), col = 'black', lwd = wid)
    #         polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1, 1, 18, 18),
    #                 col = alpha('lightgrey', 0.3), border = NA)
    #         polygon(x = c(sub.true$start[r], sub.true$end[r], sub.true$end[r], sub.true$start[r]), y = c(1-hit, 1-hit, 1+hit, 1+hit),
    #                 col = 'black', border = NA)
    #       }
    #       y <- 2
    #       for(c in c(5, 10, 15, 30, 50)){
    #         temp <- sub.pl[sub.pl$covg == c,]
    #         if(nrow(temp) > 0){
    #           for(r in 1:nrow(temp)){
    #             # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = pl.col, lwd = wid)
    #             polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
    #                     col = pl.col, border = NA)
    #           }
    #         }
    #         y <- y+1
    #       }
    #       for(c in c(5, 10, 15, 30, 50)){
    #         temp <- sub.gt[sub.gt$covg == c,]
    #         if(nrow(temp) > 0){
    #           for(r in 1:nrow(temp)){
    #             # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = gt.col, lwd = wid)
    #             polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
    #                     col = gt.col, border = NA)
    #           }
    #         }
    #         y <- y+1
    #       }
    #       for(c in c(5, 10, 15, 30, 50)){
    #         temp <- sub.pk[sub.pk$covg == c,]
    #         if(nrow(temp) > 0){
    #           for(r in 1:nrow(temp)){
    #             # lines(x = c(temp$start[r], temp$end[r]), y = c(y, y), col = plink.col, lwd = wid)
    #             polygon(x = c(temp$start[r], temp$end[r], temp$end[r], temp$start[r]), y = c(y-hit, y-hit, y+hit, y+hit),
    #                     col = plink.col, border = NA)
    #           }
    #         }
    #         y <- y+1
    #       }
    #       mtext('PLINK', side = 2, line = 3, at = 14)
    #       mtext('BCFtools\nGenotypes', side = 2, line = 3, at = 9)
    #       mtext('BCFtools\nLikelihoods', side = 2, line = 3, at = 4)
    #     }
    #     dev.off()
  }
}
  
