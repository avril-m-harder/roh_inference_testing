### Script to calculate and visualize heterozygosity in windows alongside true and called ROHs

library(scales)
library(ghibli)
library(vcfR)
`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/')

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
coverages <- cbind(c(5, 10, 15, 30, 50),
                   c('5X','10X','15X','30X','fullcovg'))

## tas dev chrom data
chroms <- read.table('data/mSarHar1.11_autosomes.txt', sep = '\t', header = TRUE)
chroms <- chroms[,c(5,9)]
colnames(chroms) <- c('chrom','length')
chroms$cum.len <- cumsum(as.numeric(chroms$length))
tot.len <- sum(chroms$length)

## Sample and chrom keys
samp.key <- read.csv('../../empirical_qc/sample_id_key.csv')
chrom.key <- read.csv('../../empirical_qc/chrom_key.csv')

## BCFtools/ROH results
gt.res <- read.csv('data/tasdev_GT_allcovg_levels.csv')
gt.res <- gt.res[gt.res$ests == 'vtrained',]
gt.res <- gt.res[gt.res$length >= 100e3,]
gt.res <- merge(gt.res, chrom.key, by.x = 'chrom', by.y = 'chrom.name')
pl.res <- read.csv('data/tasdev_PL_allcovg_levels.csv')
pl.res <- pl.res[pl.res$ests == 'vtrained',]
pl.res <- pl.res[pl.res$length >= 100e3,]
pl.res <- merge(pl.res, chrom.key, by.x = 'chrom', by.y = 'chrom.name')

## settled on default settings
plink.res <- read.csv('../../empirical_plotting_summarizing/archive_first_submission/plink_results_round1/PLINK_all_coordinates_DEFAULT_SETTINGS.csv')
plink.res$length <- plink.res$end - plink.res$start + 1

## set window size
wind.size <- 100e3

##### 2. Looping over coverage levels, read in VCF and calculate heterozygosity #####
for(cov in 1:nrow(coverages)){
  
  vcf <- read.vcfR(paste0('data/vcf_files/final_filt_tasdev_all_chroms_15samps_',coverages[cov,2],'.recode.vcf.gz'))
  gts <- extract.gt(vcf)
  chrom <- vcf@fix[,1]
  pos <- vcf@fix[,2]
  rm(vcf)
  gc()
  
  ##### 99. Plotting individual true ROHs and called ROHs for all 3 analyses and coverage levels #####
  wid <- 3 ## line widths in plot ~ or ~
  hit <- 0.1 ## polygon height
  pt.cex <- 0.4
  lwd <- 0.5
  cov.name <- coverages[cov,2]
  
  for(chrom.row in 1:nrow(chroms)){
    print(chroms$chrom[chrom.row])
    pdf(paste0('figures/',cov.name,'_',chroms$chrom[chrom.row],'_individual_true_and_called_ROHs_across_chromosome.pdf'), 
        width = 15, height = 9)
    par(mar = c(5.1, 10.1, 4.1, 2.1))
    for(i in (unique(pl.res$id))){
      id.gts <- cbind(chrom, pos, gts[,which(colnames(gts) == i)])
      id.gts <- id.gts[id.gts[,1] == chroms$chrom[chrom.row],]
      id.gts[id.gts[,3] == '0/0', 3] <- 0
      id.gts[id.gts[,3] == '0|0', 3] <- 0
      id.gts[id.gts[,3] == '1/1', 3] <- 0
      id.gts[id.gts[,3] == '1|1', 3] <- 0
      id.gts[id.gts[,3] == '0/1', 3] <- 1
      id.gts[id.gts[,3] == '1/0', 3] <- 1
      id.gts[id.gts[,3] == '0|1', 3] <- 1
      id.gts[id.gts[,3] == '1|0', 3] <- 1
      
      ## calculate window heterozygosity values
      s <- 1
      e <- s + wind.size - 1
      chrom.len <- chroms$length[chrom.row]
      OUT <- NULL
      while(e <= chrom.len){
        # print(round(e/chrom.len, digits = 2))
        mid <- e - (wind.size/2)
        het <- nrow(id.gts[which(as.numeric(id.gts[,2]) >= s &
                                 as.numeric(id.gts[,2]) <= e &
                                 as.numeric(id.gts[,3] == 1)),])
        if(is.null(het)){
          het <- 0
        }
        all <- nrow(id.gts[which(as.numeric(id.gts[,2]) >= s &
                                   as.numeric(id.gts[,2]) <= e),])
        save <- c(mid, (het/all))
        OUT <- rbind(OUT, save)
        s <- s + wind.size
        e <- s + wind.size - 1
      }
            
      sub.pl <- pl.res[pl.res$id == i & pl.res$covg == coverages[cov,1] & pl.res$chrom == chroms$chrom[chrom.row],]
      sub.gt <- gt.res[gt.res$id == i & gt.res$covg == coverages[cov,1] & gt.res$chrom == chroms$chrom[chrom.row],]
      sub.pk <- plink.res[plink.res$id == i & plink.res$covg == coverages[cov,1] & plink.res$chrom == chroms$chrom[chrom.row],]

      plot(0,0, xlim = c(1, chrom.len), ylim = c(0.85, 4.65), col = 'transparent', yaxt = 'n', xlab = 'Chromosome position (bp)', 
           ylab = '', main = '', bty = 'n')
      axis(2, at = c(1:3), labels = c('PLINK','BCFtools\nLikelihoods','BCFtools\nGenotypes'), las = 2)
      axis(2, at = c(3.5, 4, 4.5), labels = c('0','0.5','1'), las = 2)
      mtext(text = 'Heterozygosity', side = 2, las= 2, at = 4, line = 3)
      y <- 1
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
  
      OUT[,2] <- OUT[,2] + 3.5
      lines(x = c(1,chrom.len), y = c(3.5, 3.5))
      points(OUT[,1], OUT[,2], pch = 16, cex = pt.cex, col = 'springgreen4')
      # lines(OUT[,1], OUT[,2], lwd = lwd, col = 'springgreen4')
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
    #       sub.pl <- pl.res[pl.res$id == i,]
    #       sub.gt <- gt.res[gt.res$id == i,]
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
