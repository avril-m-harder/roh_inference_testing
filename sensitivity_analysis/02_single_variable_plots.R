## color settings (1 per parameter)
library(NatParksPalettes)
pal <- natparks.pals(name="DeathValley", n=7, type="discrete")
cols <- cbind(pal[c(1:3,5:7)], c('phwh','phwm','phws','phzg','phwt','phzs'))
colnames(cols) <- c('col','Variables')

# Good intro about sensitivity analysis: 
#https://complementarytraining.net/simple-sensitivity-analysis-with-r/


###### Simulated Data ######
## Input files formatted by LOCAL_03b_roh_param_analyses_preSA_plink.R script (in /simulated/R/)
setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/sensitivity_analysis/data/simulated/')

### loop over rounds/iterations (i.e., different parameter setting combinations),
### then loop over demographic scenarios within each round and analyze results
### at each level of coverage separately

## define rounds to be analyzed
rds <- list.dirs()
rds <- rds[grep('round', rds)][2]

for(rd.dir in rds){
  rd <- gsub('./', '', rd.dir)
  fns <- list.files(path = rd, pattern = 'individual')
  demos <- do.call(rbind, strsplit(fns, split = '_'))[,1]
  
  if (!file.exists(paste0('../../figures/simulated/',rd.dir))){
    dir.create(paste0('../../figures/simulated/',rd.dir))
  }
  for(d in demos){
    if (!file.exists(paste0('../../figures/simulated/',rd.dir,'/',d,'/'))){
      dir.create(paste0('../../figures/simulated/',rd.dir,'/',d,'/'))
      
      df.iter <- read.csv(paste0(rd.dir,'/',d,'_individual_froh_results_',rd,'.csv'), header = T)
      ## phzk and phzd not varied, remove from analysis
      df.iter <- df.iter[, -c(which(colnames(df.iter) %in% c('phzk','phzd','abs.diff','group')))]
      df.iter$covg <- as.numeric(gsub('x', '', df.iter$covg))
      
      g <- 1 ## set starting group number so that each coverage level has unique group numbers
      for(c in sort(unique(df.iter$covg))){
        cov.df.iter <- df.iter[df.iter$covg == c,]
        
        ## for each combination, calculate the difference between called and true f(ROH) for each individual,
        ## write results to a file, and clean up
        for(r in 1:nrow(cov.df.iter)){
          cov.df.iter$group[r] <- paste(unlist(cov.df.iter[r, c(4:10)]), collapse = '_')
        }
        key <- cbind(unique(cov.df.iter$group), seq(g, g+length(unique(cov.df.iter$group))-1))
        colnames(key) <- c('group','group.num')
        cov.df.iter <- merge(cov.df.iter, key, by = 'group')
        cov.df.iter <- cov.df.iter[,-1]
        cov.df.iter$group.num <- as.numeric(cov.df.iter$group.num)
        OUT <- NULL
        for(g in unique(cov.df.iter$group.num)){
          sub <- cov.df.iter[cov.df.iter$group.num == g,]
          m.abs <- mean(abs(sub$true.froh - sub$call.froh))
          m <- mean(sub$true.froh - sub$call.froh) ## if true is > called, on average, result will be positive
          sd.abs <- sd(abs(sub$true.froh - sub$call.froh))
          sd <- sd(sub$true.froh - sub$call.froh)
          save <- c(c, g, m, sd, m.abs, sd.abs)
          OUT <- rbind(OUT, save)
        }
        colnames(OUT) <- c('covg','group.num','mean.true.minus.call','sd','mean.abs.diff','sd.abs.diff')
        if(c == sort(unique(df.iter$covg))[1]){
          write.table(OUT, paste0(rd.dir,'/',d,'_',rd,'_true_vs_called_fROH_indiv_data.txt'),
                      row.names = FALSE, quote = FALSE, sep = '\t')
        } else{
          write.table(OUT, paste0(rd.dir,'/',d,'_',rd,'_true_vs_called_fROH_indiv_data.txt'),
                      row.names = FALSE, quote = FALSE, sep = '\t', append = TRUE)
        }
        rm(OUT, m, m.abs, r, save, sd, sd.abs, sub)
        
        # Mean and SD for each ind (calculated over all settings combos with ROH calls for each ind)
        mean.iter1 <- ddply(cov.df.iter, "id", summarise, meanTrue=mean(true.froh), sdTrue=sd(true.froh),
                            meanCall=mean(call.froh), sdCall=sd(call.froh))
        
        ## may want data like this instead of summarizing individual results?
        # summary(lm(meanCall ~ meanTrue, mean.iter1)) # Adjusted R-squared:  0.8183 ; p-value: < 2.2e-16
        
        # Standardized Rank Regression Coefficients
        inds <- unique(cov.df.iter$id)
        
        srcOriginal1 <- NULL
        srcError1 <- NULL
        srcBias1 <- NULL
        for (ind in inds)
        {
          df <- cov.df.iter[which(cov.df.iter$id ==ind),]
          if(length(unique(df$group.num)) > 8){
            
            targets <- df[,3]
            ## detect predictors for each round based on which are variable
            CS <- NULL
            for(p in c(5:10)){
              if(length(unique(df[,p])) > 1){
                CS <- c(CS, p)
              }
            }
            predictors <- df[,c(CS)]
            srcInd <- sensitivity::src(predictors, targets, rank = T, logistic = TRUE, nboot = 100, conf = 0.95)
            df3 <- c(ind,as.array(t(srcInd$SRC))[1,])
            df4 <- c(ind,as.array(t(srcInd$SRC))[3,])
            df5 <- c(ind,as.array(t(srcInd$SRC))[2,])
            srcOriginal1 <- rbind(srcOriginal1,df3)
            srcError1 <- rbind(srcError1,df4)
            srcBias1 <- rbind(srcBias1,df5)
          } else{
            next
          }
        }
        
        srcOriginal1 <- as.data.frame(srcOriginal1)
        srcError1 <- as.data.frame(srcError1)
        srcBias1 <- as.data.frame(srcBias1)
        
        colnames(srcOriginal1)[1] <- "id"
        colnames(srcBias1)[1] <- "id"
        print(paste0(rd,' - ',c,'X: ',length(unique(srcOriginal1$id)),' samples'))
        
        # Plot SRCs for each variable
        origMelt1 <- melt(data = srcOriginal1, id.vars="id",
                          measure.vars = colnames(srcOriginal1)[colnames(srcOriginal1) %notin% c('id')],
                          variable.name = "Variables",
                          value.name = "SRC")
        
        ### base R figure
        origMelt1 <- merge(origMelt1, cols, by = 'Variables')
        max <- length(unique(origMelt1$Variables))
        text.size <- 1.25
        
        
        shrink <- 400 ## higher # here ==> narrower x-direction spread for points
        alph <- 0.4
        pt.size <- 0.8
        xmin <- 0.85
        xmax <- max + 0.15
        
        if(length(unique(origMelt1$Variables)) > 2){
          pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_paramvsSRC.pdf'), width = 6, height = 5)
        } else{
          pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_paramvsSRC.pdf'), width = 3, height = 5)
          xmin <- 0.65
          xmax <- max + 0.35
        }
        x <- 1
        plot(0,0, xlim = c(xmin, xmax), ylim = c(-1, 1), xaxt = 'n', xlab = 'Parameter', ylab = 'Standardized regression coefficient', col = 'transparent', cex.axis = text.size, cex.lab = text.size, main = paste0(rd,' - ',c,'X'))
        axis(1, at = c(1:max), labels = unique(origMelt1$Variables), cex.axis = text.size)
        abline(h = 0, lty = 2)
        for(v in unique(origMelt1$Variables)){
          sub <- origMelt1[origMelt1$Variables == v,]
          f <- sample(c(-100:100), nrow(sub))
          f <- f/shrink+x
          points(f, sub$SRC, col = alpha(sub$col[1], alph), pch = 19, cex = pt.size)
          arrows(x0 = x, x1 = x, y0 = (mean(sub$SRC, na.rm = TRUE) - sd(sub$SRC, na.rm = TRUE)),
                 y1 = (mean(sub$SRC, na.rm = TRUE) + sd(sub$SRC, na.rm = TRUE)),
                 lwd = 2, col = 'black', code=3, angle=90, length=0)
          points(x, mean(sub$SRC), pch = 23, col = 'black', bg = sub$col[1], cex = pt.size + 0.5, lwd = 2)
          
          x <- x+1
        }
        dev.off()
        
        # SRCs vs true value
        rohSRC1 <- inner_join(srcOriginal1,mean.iter1,by="id")
        
        # Plot each variable vs true FROH
        rohMelt1 <- melt(data = rohSRC1, id.vars=c("id","meanTrue"),
                         measure.vars = colnames(rohSRC1)[colnames(rohSRC1) %notin% c('id','meanTrue','sdTrue','meanCall','sdCall')],
                         variable.name = "Variables",
                         value.name = "SRC")
        
        rohMelt1 <- merge(rohMelt1, cols, by = 'Variables')
        alph <- 0.8
        pt.size <- 0.8
        txt.size <- 1.25
        
        pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_truefROH_vs_SRC.pdf'), width = 7, height = 5)
        par(mai = c(1.02,0.82,0.82,1.42))
        plot(rohMelt1$meanTrue, rohMelt1$SRC, ylim = c(-1, 1), col = 'transparent', pch = 19, xlab = substitute(paste('Mean true ',italic('F')[ROH])), ylab = 'Standardized regression coefficient', cex.axis = txt.size, cex.lab = txt.size,
             main = paste0(rd,' - ',c,'X'))
        abline(h = 0, lty = 2)
        points(rohMelt1$meanTrue, rohMelt1$SRC, col = alpha(rohMelt1$col, alph), pch = 19, cex = pt.size)
        for(v in unique(rohMelt1$Variables)){
          sub <- rohMelt1[rohMelt1$Variables == v,]
          vals <- loess.smooth(sub$meanTrue, sub$SRC, span = 0.75,
                               family = c("gaussian"), col = sub$col[1])
          lines(vals$x, vals$y, col = sub$col[1], lwd = 3)
        }
        par(xpd = TRUE)
        legend('right', lwd = 3, col = unique(rohMelt1$col), legend = unique(rohMelt1$Variables), bty = 'n', inset = -0.30, cex = txt.size)
        dev.off()
        
        
        # Plot each variable vs called FROH
        rohMelt1 <- melt(data = rohSRC1, id.vars=c("id","meanCall"),
                         measure.vars = colnames(rohSRC1)[colnames(rohSRC1) %notin% c('id','meanTrue','sdTrue','meanCall','sdCall')],
                         variable.name = "Variables",
                         value.name = "SRC")
        
        rohMelt1 <- merge(rohMelt1, cols, by = 'Variables')
        alph <- 0.8
        pt.size <- 0.8
        txt.size <- 1.25
        
        pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_callfROH_vs_SRC.pdf'), width = 7, height = 5)
        par(mai = c(1.02,0.82,0.82,1.42))
        plot(rohMelt1$meanCall, rohMelt1$SRC, ylim = c(-1, 1), col = 'transparent', pch = 19, xlab = substitute(paste('Mean called ',italic('F')[ROH])), ylab = 'Standardized regression coefficient', cex.axis = txt.size, cex.lab = txt.size,
             main = paste0(rd,' - ',c,'X'))
        abline(h = 0, lty = 2)
        points(rohMelt1$meanCall, rohMelt1$SRC, col = alpha(rohMelt1$col, alph), pch = 19, cex = pt.size)
        for(v in unique(rohMelt1$Variables)){
          sub <- rohMelt1[rohMelt1$Variables == v,]
          vals <- loess.smooth(sub$meanCall, sub$SRC, span = 0.75,
                               family = c("gaussian"), col = sub$col[1])
          lines(vals$x, vals$y, col = sub$col[1], lwd = 3)
        }
        par(xpd = TRUE)
        legend('right', lwd = 3, col = unique(rohMelt1$col), legend = unique(rohMelt1$Variables), bty = 'n', inset = -0.30, cex = txt.size)
        dev.off()
        
        # Plot true vs. called f(ROH)
        pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_true_v_called_fROH.pdf'), width = 7.05, height = 5)
        par(mai = c(1.02,0.87,0.82,1.02))
        plot(rohSRC1$meanTrue, rohSRC1$meanCall, col = 'transparent', pch = 19, xlab = substitute(paste('True ',italic('F')[ROH])),
             ylab = substitute(paste('Mean called ',italic('F')[ROH])), cex.axis = txt.size, cex.lab = txt.size,
             main = paste0(rd,' - ',c,'X'), ylim = c(min(rohSRC1$meanCall - rohSRC1$sdCall), max(rohSRC1$meanCall + rohSRC1$sdCall)))
        # plot(rohSRC1$meanTrue, rohSRC1$meanCall, col = 'transparent', pch = 19, xlab = substitute(paste('True ',italic('F')[ROH])),
        #      ylab = substitute(paste('Mean called ',italic('F')[ROH])), cex.axis = txt.size, cex.lab = txt.size,
        #      main = paste0(rd,' - ',c,'X'), ylim = c(0, 0.35), yaxt = 'n')
        #   axis(2, at = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35), labels = c('0.0','','0.1','','0.2','','0.3',''), cex.lab = txt.size, cex.axis = txt.size)
        
        abline(0, 1, lty = 2)
        points(rohSRC1$meanTrue, rohSRC1$meanCall, col = alpha(pal[4], 1), pch = 19)
        arrows(x0 = rohSRC1$meanTrue, x1 = rohSRC1$meanTrue, y0 = rohSRC1$meanCall - rohSRC1$sdCall, y1 = rohSRC1$meanCall + rohSRC1$sdCall,
               lwd = 1, col = pal[4], code=3, angle=90, length=0)
        dev.off()
      }
    }
    if(file.exists(paste0('../../figures/simulated/',rd.dir,'/',d,'/'))){
      print(paste0('Already processed: ',rd,' - ',d))
    }
  }
}