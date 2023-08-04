library(scales)

setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/empirical/data/')

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

bin.names <- c('short ROHs','intermediate ROHs','long ROHs','very long ROHs')

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
colnames(dat) <- c('method','id','covg','ests','froh','froh.1','froh.2','froh.3','froh.4')

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
              col = alpha('black', 1), pch = 19, cex = 0.75)
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
                 col = alpha('black', 1), pch = 19, cex = 0.75)
        }
        x <- x+2
      }
      n <- n+1
    }
    dev.off()
}




