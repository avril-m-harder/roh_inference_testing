library(ggplot2)
library(plyr)
library(dplyr)
library(DataCombine)
library(sensitivity)
library(reshape2)
library(TeachingDemos)
library(scales)
`%notin%` <- Negate(`%in%`)

## color settings (1 per parameter)
library(NatParksPalettes)
pal <- natparks.pals(name="DeathValley", n=7, type="discrete")
cols <- cbind(pal[c(1:3,5:7)], c('phwh','phwm','phws','phzg','phwt','phzs'))
colnames(cols) <- c('col','Variables')


###### Simulated Data ######
## Input files formatted by LOCAL_03b_roh_param_analyses_preSA_plink.R script (in /simulated/R/)
setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/sensitivity_analysis/data/simulated/')

## define rounds to be analyzed
rds <- list.dirs()
rds <- rds[grep('round', rds)][1]

## define PLINK default settings for parameters
defaults <- as.data.frame(cbind(c('phwh','phwm','phws','phzg','phwt','phzs'),
                  c(1, 5, 50, 1000, 0.05, 100)))
colnames(defaults) <- c('param','value')

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
    
    for(V in VAR){
      var.df.iter <- df.iter
      vars.to.set <- defaults[defaults$param != V, 1]
      vals.to.set <- defaults[defaults$param != V, 2]
      for(v in 1:length(vars.to.set)){
        colm <- which(colnames(df.iter) == vars.to.set[v])
        var.df.iter <- var.df.iter[var.df.iter[, colm] == vals.to.set[v],]
      }
      
      for(c in unique(var.df.iter$covg)){
        cov.df.iter <- var.df.iter[var.df.iter$covg == c,]
        
        # Mean and SD for each ind (calculated over all settings combos with ROH calls for each ind)
        mean.iter1 <- ddply(cov.df.iter, "id", summarise, meanTrue=mean(true.froh), sdTrue=sd(true.froh),
                            meanCall=mean(call.froh), sdCall=sd(call.froh))
        
        ## True vs. called f(ROH)
        pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_',V,'_varsettingcomps_truevscalledfroh.pdf'), width = 6, height = 5)
        par(mar = c(5.1, 5.1, 4.1, 6.1), xpd = FALSE)
        
        x <- 1
        colour <- cols[cols[,2] == V, 1]
        pchs <- c(16, 17, 15, 18, 3, 4)
        sets <- sort(unique(cov.df.iter[,V]))
        plot.min <- min(cov.df.iter$true.froh, cov.df.iter$call.froh)
        plot.max <- max(cov.df.iter$true.froh, cov.df.iter$call.froh)
        max <- length(unique(cov.df.iter[,V]))
        text.size <- 1.25
        alph <- 0.6
        pt.size <- 1
        
        plot(0,0, xlim = c(plot.min, plot.max), ylim = c(plot.min, plot.max), 
             xlab = substitute(paste('True ',italic('F')[ROH])), 
             ylab = substitute(paste('Called ',italic('F')[ROH])), col = 'transparent', 
             cex.axis = text.size, cex.lab = text.size, main = paste0(rd,' - ',c,'X'))
          abline(0, 1, lty = 2)
          for(v in 1:max){
            sub <- cov.df.iter[which(cov.df.iter[,V] == sets[v]),]
            points(sub$true.froh, sub$call.froh, col = alpha(colour, alph), pch = pchs[v], cex = pt.size)
            if(nrow(sub) > 20){
              suppressWarnings(clipplot(abline(lm(sub$call.froh ~ sub$true.froh), col = colour, lty = v), 
                       xlim = c(min(sub$true.froh), max(sub$true.froh))))
            }
          }
          par(xpd = TRUE)
          legend('right', inset = c(-0.26, 0), lty = c(1:max), pch = pchs[c(1:max)], 
                 legend = sets, col = colour, title = paste0(V,'\ndefault = ',defaults[defaults$param == V, 2]), bty = 'n')
        dev.off()
  
        ## Variation in called f(ROH) across variable settings
        x.max <- max + 0.25
        pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_',V,'_varsettingcomps_callfROH_by_value.pdf'), width = 6, height = 5)
        par(mar = c(5.1, 5.1, 4.1, 2.1), xpd = FALSE)
  
        plot(0,0, xlim = c(0.75, x.max), ylim = c(plot.min, plot.max), col = 'transparent', 
             cex.axis = text.size, cex.lab = text.size, main = paste0(rd,' - ',c,'X\n',V), xlab = '', xaxt = 'n',
             ylab = substitute(paste('Called ',italic('F')[ROH])))
          axis(1, at = c(1:max), cex.axis = text.size, labels = sets)
          for(v in 1:max){
            sub <- cov.df.iter[which(cov.df.iter[,V] == sets[v]),]
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
       col <- which(colnames(cov.df.iter) == V)
       for(p in sets){
         OUT <- suppressWarnings(cbind(OUT, cov.df.iter[which(cov.df.iter[,col] == p), 'call.froh']))
       }
       colnames(OUT) <- c('id', sets)
       pdf(paste0('../../figures/simulated/',rd.dir,'/',d,'/',d,'_simulated_',rd,'_',c,'X_',V,'_varsettingcomps_setting_correlations.pdf'), width = 6, height = 5)
       pairs(OUT[,which(colnames(OUT) != 'id')], pch = 16, col = alpha(colour, alph), 
             upper.panel = panel.cor,
             gap=0, row1attop = FALSE)
       dev.off()
      }
    }
  }  
  if(file.exists(paste0('../../figures/simulated/',rd.dir,'/',d,'/'))){
    print(paste0('Already processed: ',rd,' - ',d))
  }
}
