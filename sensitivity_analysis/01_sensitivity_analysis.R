###########################################################################
###                          Samarth Mathur, PhD                        ###
###                       The Ohio State University                     ###
###                                                                     ###
###     Date Created: 07/11/22                  Last Modified: 07/11/22 ###
###########################################################################
###########################################################################
###                   ROH_sensitivity.R  		                            ###
###########################################################################

#### PREREQUISITES #####
# load packages
library(ggplot2)
library(plyr)
library(dplyr)
library(DataCombine)
library(sensitivity)
library(reshape2)
library(TeachingDemos)
`%notin%` <- Negate(`%in%`)

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
    }
      
    df.iter <- read.csv(paste0(rd.dir,'/',d,'_individual_froh_results_',rd,'.csv'), header = T)
    ## phzk and phzd not varied, remove from analysis
    df.iter <- df.iter[, -c(which(colnames(df.iter) %in% c('phzk','phzd','abs.diff','group')))]
    df.iter$covg <- as.numeric(gsub('x', '', df.iter$covg))
    
    ## check how many variables are.. variable
    ct <- 0
    VAR <- NULL
    for(s in c('phwh','phwm','phws','phzg','phwt','phzs')){
      if(length(table(df.iter[,s])) > 1){
        ct <- ct+1
        VAR <- c(VAR, s)
      }
    }
    ## if >1 parameter had variable settings, go ahead and do the sensitivity analysis.
    if(ct > 1){
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
        
        if(length(unique(origMelt1$Variables)) > 3){
          pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_paramvsSRC.pdf'), width = 6, height = 5)
        } else if(length(unique(origMelt1$Variables)) == 3){
          pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_paramvsSRC.pdf'), width = 4.5, height = 5)
          xmin <- 0.65
          xmax <- max + 0.35
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
    ## otherwise, just make some plots to examine the single variable parameter setting
    } else{
      for(c in unique(df.iter$covg)){
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
      
      ## True vs. called f(ROH)
      pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_truevscalledfroh.pdf'), width = 6, height = 5)
      par(mar = c(5.1, 5.1, 4.1, 6.1), xpd = FALSE)
      
      x <- 1
      colour <- cols[cols[,2] == s, 1]
      pchs <- c(16, 17, 15, 18, 3, 4)
      sets <- sort(unique(cov.df.iter[,s]))
      plot.min <- min(cov.df.iter$true.froh, cov.df.iter$call.froh)
      plot.max <- max(cov.df.iter$true.froh, cov.df.iter$call.froh)
      max <- length(unique(cov.df.iter[,s]))
      text.size <- 1.25
      alph <- 0.6
      pt.size <- 1
      
      plot(0,0, xlim = c(plot.min, plot.max), ylim = c(plot.min, plot.max), 
           xlab = substitute(paste('True ',italic('F')[ROH])), 
           ylab = substitute(paste('Called ',italic('F')[ROH])), col = 'transparent', 
           cex.axis = text.size, cex.lab = text.size, main = paste0(rd,' - ',c,'X'))
        abline(0, 1, lty = 2)
        for(v in 1:max){
          sub <- cov.df.iter[which(cov.df.iter[,s] == sets[v]),]
          points(sub$true.froh, sub$call.froh, col = alpha(colour, alph), pch = pchs[v], cex = pt.size)
          suppressWarnings(clipplot(abline(lm(sub$call.froh ~ sub$true.froh), col = colour, lty = v), 
                   xlim = c(min(sub$true.froh), max(sub$true.froh))))
        }
        par(xpd = TRUE)
        legend('right', inset = c(-0.26, 0), lty = c(1:max), pch = pchs[c(1:max)], 
               legend = unique(cov.df.iter[,s]), col = colour, title = s, bty = 'n')
      dev.off()

      ## Variation in called f(ROH) across variable settings
      x.max <- max + 0.25
      pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_callfROH_by_value.pdf'), width = 6, height = 5)
      par(mar = c(5.1, 5.1, 4.1, 2.1), xpd = FALSE)

      plot(0,0, xlim = c(0.75, x.max), ylim = c(plot.min, plot.max), col = 'transparent', 
           cex.axis = text.size, cex.lab = text.size, main = paste0(rd,' - ',c,'X\n',s), xlab = '', xaxt = 'n',
           ylab = substitute(paste('Called ',italic('F')[ROH])))
        axis(1, at = c(1:max), cex.axis = text.size, labels = sets)
        for(v in 1:max){
          sub <- cov.df.iter[which(cov.df.iter[,s] == sets[v]),]
          points(jitter(rep(v, nrow(sub)), amount = 0.15), sub$call.froh, col = alpha(colour, alph), pch = pchs[v], cex = pt.size)
          points(v, mean(sub$call.froh), pch = 23, col = 'black', bg = colour, cex = pt.size + 0.5, lwd = 2)
        }
      dev.off()

     ## Correlation of called f(ROH) values across parameter values
     panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
     {
       usr <- par("usr"); on.exit(par(usr))
       par(usr = c(0, 1, 0, 1))
       r <- abs(cor(x, y))
       txt <- format(c(r, 0.123456789), digits = digits)[1]
       txt <- paste0(prefix, txt)
       if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
       text(0.5, 0.5, txt, cex = 1.5)
     }
     
     OUT <- NULL
     cov.df.iter <- cov.df.iter[order(cov.df.iter$id),]
     OUT <- cbind(OUT, unique(cov.df.iter$id))
     col <- which(colnames(cov.df.iter) == s)
     for(p in sets){
       OUT <- cbind(OUT, cov.df.iter[which(cov.df.iter[,col] == p), 'call.froh'])
     }
     colnames(OUT) <- c('id', sets)
     pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_setting_correlations.pdf'), width = 6, height = 5)
     pairs(OUT[,which(colnames(OUT) != 'id')], pch = 16, col = alpha(colour, alph), 
           upper.panel = panel.cor,
           gap=0, row1attop = FALSE)
     dev.off()
    }
  }  
  if(file.exists(paste0('../../figures/simulated/',rd.dir,'/',d,'/'))){
    print(paste0('Already processed: ',rd,' - ',d))
  }
  }
}

##### Empirical Data (2023: needs to be revised!) #####
# setwd('/Users/Avril/Desktop/roh_param_project/iterative_sa_steps/empirical/')
# 
# ### loop over rounds/iterations (i.e., different parameter setting combinations) and analyze results 
# ### at each level of coverage separately
# 
# ## define rounds to be analyzed
# # rds <- c('round1','round2a','round2b','round2c','round2d','round2e')
# rds <- c('round1','round2a','round2c','round2d','round2e')
# 
# for(rd in rds){
#   df.iter <- read.csv(paste0(rd,'/individual_froh_results_',rd,'.csv'), header = T)
#   ## phzk and phzd not varied, remove from analysis
#   df.iter <- df.iter[, -c(which(colnames(df.iter) %in% c('phzk','phzd','abs.diff','group')))]
#   
#   g <- 1 ## set starting group number so that each coverage level has unique group numbers
#   for(c in sort(unique(df.iter$covg))){
#     cov.df.iter <- df.iter[df.iter$covg == c,]
#     
#     ## for each combination, calculate the difference between called and true f(ROH) for each individual,
#     ## write results to a file, and clean up
#     for(r in 1:nrow(cov.df.iter)){
#       cov.df.iter$group[r] <- paste(unlist(cov.df.iter[r, c(3:9)]), collapse = '_')
#     }
#     key <- cbind(unique(cov.df.iter$group), seq(g, g+length(unique(cov.df.iter$group))-1))
#     colnames(key) <- c('group','group.num')
#     cov.df.iter <- merge(cov.df.iter, key, by = 'group')
#     cov.df.iter <- cov.df.iter[,-1]
#     cov.df.iter$group.num <- as.numeric(cov.df.iter$group.num)
#     
#     # Mean and SD for each ind (calculated over all settings combos with ROH calls for each ind)
#     mean.iter1 <- ddply(cov.df.iter, "id", summarise, meanCall=mean(call.froh), sdCall=sd(call.froh))
#     
#     ## may want data like this instead of summarizing individual results?
#     # summary(lm(meanCall ~ meanTrue, mean.iter1)) # Adjusted R-squared:  0.8183 ; p-value: < 2.2e-16
#     
#     # Standardized Rank Regression Coefficients
#     inds <- unique(cov.df.iter$id)
#     
#     srcOriginal1 <- NULL
#     srcError1 <- NULL
#     srcBias1 <- NULL
#     for (ind in inds)
#     {
#       df <- cov.df.iter[which(cov.df.iter$id ==ind),]
#       if(length(unique(df$group.num)) > 5){
#         
#         targets <- df[,2]
#         ## detect predictors for each round based on which are variable
#         CS <- NULL
#         for(p in c(4:9)){
#           if(length(unique(df[,p])) > 1){
#             CS <- c(CS, p)
#           }
#         }
#         predictors <- df[,c(CS)]
#         srcInd <- sensitivity::src(predictors, targets, rank = T, logistic = TRUE, nboot = 100, conf = 0.95)
#         df3 <- c(ind,as.array(t(srcInd$SRC))[1,])
#         df4 <- c(ind,as.array(t(srcInd$SRC))[3,])
#         df5 <- c(ind,as.array(t(srcInd$SRC))[2,])
#         srcOriginal1 <- rbind(srcOriginal1,df3)
#         srcError1 <- rbind(srcError1,df4)
#         srcBias1 <- rbind(srcBias1,df5)
#       } else{
#         next
#       }
#       
#     }
#     
#     srcOriginal1 <- as.data.frame(srcOriginal1)
#     srcError1 <- as.data.frame(srcError1)
#     srcBias1 <- as.data.frame(srcBias1)
#     
#     colnames(srcOriginal1)[1] <- "id"
#     colnames(srcBias1)[1] <- "id"
#     print(paste0(rd,' - ',c,'X: ',length(unique(srcOriginal1$id)),' samples'))
#     
#     # Plot SRCs for each variable
#     origMelt1 <- melt(data = srcOriginal1, id.vars="id",
#                       measure.vars = colnames(srcOriginal1)[colnames(srcOriginal1) %notin% c('id')],
#                       variable.name = "Variables",
#                       value.name = "SRC")
#     
#     ### base R figure
#     origMelt1 <- merge(origMelt1, cols, by = 'Variables')
#     origMelt1$SRC <- as.numeric(origMelt1$SRC)
#     max <- length(unique(origMelt1$Variables))
#     text.size <- 1.25
#     shrink <- 400 ## higher # here ==> narrower x-direction spread for points
#     alph <- 0.4
#     pt.size <- 0.8
#     xmin <- 0.85
#     xmax <- max + 0.15
#     
#     if(length(unique(origMelt1$Variables)) > 2){
#       pdf(paste0('figures/empirical_',rd,'_',c,'X_paramvsSRC.pdf'), width = 6.25, height = 5)
#     } else{
#       pdf(paste0('figures/empirical_',rd,'_',c,'X_paramvsSRC.pdf'), width = 3, height = 5)
#       xmin <- 0.65
#       xmax <- max + 0.35
#     }
#     x <- 1
#     plot(0,0, xlim = c(xmin, xmax), ylim = c(-1, 1), xaxt = 'n', xlab = 'Parameter', ylab = 'Standardized regression coefficient', col = 'transparent', cex.axis = text.size, cex.lab = text.size, main = paste0(rd,' - ',c,'X'))
#       axis(1, at = c(1:max), labels = unique(origMelt1$Variables), cex.axis = text.size)
#       abline(h = 0, lty = 2)
#       for(v in unique(origMelt1$Variables)){
#         sub <- origMelt1[origMelt1$Variables == v,]
#         f <- sample(c(-100:100), nrow(sub))
#         f <- f/shrink+x
#         points(f, sub$SRC, col = alpha(sub$col[1], alph), pch = 19, cex = pt.size)
#         arrows(x0 = x, x1 = x, y0 = (mean(sub$SRC, na.rm = TRUE) - sd(sub$SRC, na.rm = TRUE)),
#                y1 = (mean(sub$SRC, na.rm = TRUE) + sd(sub$SRC, na.rm = TRUE)),
#                lwd = 2, col = 'black', code=3, angle=90, length=0)
#         points(x, mean(sub$SRC), pch = 23, col = 'black', bg = sub$col[1], cex = pt.size + 0.5, lwd = 2)
#         
#         x <- x+1
#       }
#     dev.off()
#     
#     # SRCs vs true value
#     rohSRC1 <- inner_join(srcOriginal1,mean.iter1,by="id")
#     
#     # Plot each variable vs true FROH
#     rohMelt1 <- melt(data = rohSRC1, id.vars=c("id","meanCall"),
#                      measure.vars = colnames(rohSRC1)[colnames(rohSRC1) %notin% c('id','meanCall','sdCall')],
#                      variable.name = "Variables",
#                      value.name = "SRC")
#     
#     rohMelt1 <- merge(rohMelt1, cols, by = 'Variables')
#     alph <- 0.8
#     pt.size <- 0.8
#     txt.size <- 1.25
#     
#     pdf(paste0('figures/empirical_',rd,'_',c,'X_callfROH_vs_SRC.pdf'), width = 7, height = 5)
#     par(mai = c(1.02,0.82,0.82,1.42))
#     plot(rohMelt1$meanCall, rohMelt1$SRC, ylim = c(-1, 1), col = 'transparent', pch = 19, xlab = substitute(paste('Mean called ',italic('F')[ROH])), ylab = 'Standardized regression coefficient', cex.axis = txt.size, cex.lab = txt.size,
#          main = paste0(rd,' - ',c,'X'))
#       abline(h = 0, lty = 2)
#       points(rohMelt1$meanCall, rohMelt1$SRC, col = alpha(rohMelt1$col, alph), pch = 19, cex = pt.size)
#       for(v in unique(rohMelt1$Variables)){
#         sub <- rohMelt1[rohMelt1$Variables == v,]
#         vals <- loess.smooth(sub$meanCall, sub$SRC, span = 0.75,
#                              family = c("gaussian"), col = sub$col[1])
#         lines(vals$x, vals$y, col = sub$col[1], lwd = 3)
#       }
#       par(xpd = TRUE)
#       legend('right', lwd = 3, col = unique(rohMelt1$col), legend = unique(rohMelt1$Variables), bty = 'n', inset = -0.30, cex = txt.size)
#     dev.off()
#   }
# }
# 
