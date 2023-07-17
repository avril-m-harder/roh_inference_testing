library(scales)
library(grDevices)
`%notin%` <- Negate(`%in%`)

setwd('/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/plink_results_round1/')

##### Read in known/true heterozygosity + ROH information #####
true.rohs <- read.table('true_indiv_data/true_roh_coords.txt')
colnames(true.rohs) <- c('id','start','end','length')
true.rohs <- true.rohs[true.rohs$length >= 100000,] ### only keeping true ROHs >= 100 kb because that's all we're evaluating the ability to call
## create unique ROH ID for linking true ROHs to called ROHs
ns <- c(1:nrow(true.rohs))
true.rohs$true.roh.id <- ns
chrom.len <- 30e6

##### Read in ROH calling results #####
## PLINK (written from 02b_summarize_plink_output.R)
plink.res <- read.table('PLINK_all_coordinates.txt')
colnames(plink.res) <- c('id','start','end','n.snps','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','pop.size','covg')
plink.res$called.roh.id <- c(1:nrow(plink.res))

## define parameter settings
phwh=c(0, 1, 2)           # Values for -homozyg-window-het
phwm=c(2, 5, 50)          # Values for -homozyg-window-missing
phws=c(50, 100, 1000)     # Values for -homozyg-window-snp
phzd=c(50)                # Values for -homozyg-density
phzg=c(500, 1000)         # Values for -homozyg-gap
phwt=c(0.01, 0.05, 0.1)   # Values for -homozyg-window-threshold
phzs=c(10, 100, 1000)     # Values for -homozyg-snp
phzk=c(100)               # Values for -homozyg-kb
pop.size <- c(100)            # Values for sample size
covg <- c(5,10,15,30,50)              # Values for mean coverage


##### 2. ROH identification -- PLINK results #####
## Output: for each true ROH, identify all overlapping called ROHs. Each instance
## of overlap will occupy one line in the output matrix. The output matrix will contain:
## roh.id for true ROH, id/pop.size/covg/start/end for called ROH, and length of overlap.

## For each true ROH, calculate...
## the proportion overlapping a called ROH
# PLINK.OUT <- NULL
# for(i in unique(plink.res$id)){              ## for each individual,
#   print(i)
#   true.sub <- true.rohs[true.rohs$id == i,]  ## subset true ROH data,
#   sub.plink <- plink.res[plink.res$id == i,]  ## subset PLINK results
#   for(r in 1:nrow(true.sub)){                     ## loop over true ROHs,
#     s <- true.sub$start[r]                        ## save true start,
#     e <- true.sub$end[r]                          ## and true end.
#     ## check for overlaps in PLINK results
#     ## called ROHs beginning outside of true ROH, ending inside
#     if(nrow(sub.plink[sub.plink$start < s & sub.plink$end >= s & sub.plink$end <= e,]) > 0){
#       temp <- sub.plink[sub.plink$start < s & sub.plink$end >= s & sub.plink$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         true.len <- temp$end[t] - s + 1                                  ## length of overlap with the focal true ROH
#         save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t],
#                              temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t],
#                              len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t]))
#         PLINK.OUT <- rbind(PLINK.OUT, save)
#       }
#     }
#     write.table(PLINK.OUT, 'PLINK_overlap_results.txt',
#                 sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
#     PLINK.OUT <- NULL
#     ## called ROHs beginning inside of a true ROH, ending outside
#     if(nrow(sub.plink[sub.plink$start >= s & sub.plink$start <= e & sub.plink$end > e,]) > 0){
#       temp <- sub.plink[sub.plink$start >= s & sub.plink$start <= e & sub.plink$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         true.len <- e - temp$start[t] + 1                                ## length of overlap with the focal true ROH
#         save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t],
#                              temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t],
#                              len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t]))
#         PLINK.OUT <- rbind(PLINK.OUT, save)
#       }
#     }
#     write.table(PLINK.OUT, 'PLINK_overlap_results.txt',
#                 sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
#     PLINK.OUT <- NULL
#     ## called ROHs completely covering a true ROH
#     if(nrow(sub.plink[sub.plink$start < s & sub.plink$end > e,]) > 0){
#       temp <- sub.plink[sub.plink$start < s & sub.plink$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         true.len <- e - s + 1                                            ## length of overlap with the focal true ROH,
#         save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t],
#                              temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t],
#                              len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t]))
#         PLINK.OUT <- rbind(PLINK.OUT, save)
#       }
#     }
#     write.table(PLINK.OUT, 'PLINK_overlap_results.txt',
#                 sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
#     PLINK.OUT <- NULL
#     ## called ROHs completely within a true ROH
#     if(nrow(sub.plink[sub.plink$start >= s & sub.plink$end <= e,]) > 0){
#       temp <- sub.plink[sub.plink$start >= s & sub.plink$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each called ROH,
#         len <- temp$end[t] - temp$start[t] + 1                           ## calculate its length,
#         true.len <- len                                                  ## length of overlap with the focal true ROH,
#         save <- as.numeric(c(temp$id[t], temp$pop.size[t], temp$covg[t], temp$phwh[t], temp$phwm[t], temp$phws[t],
#                              temp$phzd[t], temp$phzg[t], temp$phwt[t], temp$phzs[t], temp$phzk[t],
#                              len, true.len, true.sub$true.roh.id[r], temp$called.roh.id[t], temp$n.snps[t]))
#         PLINK.OUT <- rbind(PLINK.OUT, save)
#       }
#     }
#     write.table(PLINK.OUT, 'PLINK_overlap_results.txt',
#                 sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
#     PLINK.OUT <- NULL
#   }
# }
# rm(PLINK.OUT)
plink.out <- read.table('PLINK_overlap_results.txt',
                        sep = '\t')
colnames(plink.out) <- c('id','pop.size','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','len','true.len','true.roh.id','called.roh.id','n.snps')
plink.out <- plink.out[plink.out$len >= 100000,]


##### 3. Plotting results #####
##### >>> Plot absolute difference between true and called f(ROH) #####
# OUT <- NULL
# write.table(OUT, 'individual_froh_results.txt',
#             sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
# ct <- 1
# res <- 0
# # for(pop in pop.size){
# for(pop in 100){
#   for(cov in covg){
#     for(a in phwh){
#       for(b in phwm){
#         for(c in phws){
#           for(d in phzd){
#             for(e in phzg){
#               for(f in phwt){
#                 for(g in phzs){
#                   for(h in phzk){
#                     print(ct/2430)
#                     sub <- plink.out[plink.out$pop.size == pop & plink.out$covg == cov & plink.out$phwh == a &
#                                      plink.out$phwm == b & plink.out$phws == c & plink.out$phzd == d &
#                                      plink.out$phzg == e & plink.out$phwt == f & plink.out$phzs == g &
#                                      plink.out$phzk == h,]
#                     if(nrow(sub) > 0){
#                       res <- res+1
#                       for(i in unique(sub$id)){
#                         temp <- sub[sub$id == i, c('len','called.roh.id')]
#                         temp <- temp[!duplicated(temp),]
#                         true.froh <- sum(true.rohs[true.rohs$id == i, 'length'])/30e6
#                         call.froh <- sum(temp$len)/30e6
#                         save <- c(i, pop, cov, a, b, c, d, e, f, g, h, true.froh, call.froh)
#                         OUT <- rbind(OUT, save)
#                       }
#                     }
#                     ct <- ct+1
#                     write.table(OUT, 'individual_froh_results.txt',
#                                 sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
#                     OUT <- NULL
#                   }
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }
# OUT <- read.table('individual_froh_results.txt')
# colnames(OUT) <- c('id','pop.size','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','true.froh','call.froh')
# froh.stats <- as.data.frame(OUT)
# froh.stats$abs.diff <- abs(froh.stats$call.froh - froh.stats$true.froh)
# froh.stats$group <- paste0(froh.stats$pop.size,'-',froh.stats$covg,'-',froh.stats$phwh,'-',froh.stats$phwm,'-',froh.stats$phws,'-',froh.stats$phzd,'-',froh.stats$phzg,'-',froh.stats$phwt,'-',froh.stats$phzs,'-',froh.stats$phzk)
# froh.stats <- froh.stats[froh.stats$pop.size == 100,]
# froh.stats <- froh.stats[,c('id','true.froh','call.froh','abs.diff','group','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk')]
# write.csv(froh.stats, 'individual_froh_results.csv', row.names = FALSE)
froh.stats <- read.csv('individual_froh_results.csv')

## add f(ROH) difference with sign (i.e., directionality)
froh.stats$froh.diff <- froh.stats$call.froh - froh.stats$true.froh
#### if positive, overestimated
#### if negative, underestimated


##### Writing files for sensitivity analyses (Samarth) #####
## calculate population-level data for sensitivity analyses;
## only considering combinations of settings for which there are ROHs in all 100 samples
OUT1 <- NULL
ct <- 1
for(g in unique(froh.stats$group)){
  print(ct/length(unique(froh.stats$group)))
  sub <- froh.stats[froh.stats$group == g,]
  if(length(unique(sub$id)) == 100){
    mean.true <- mean(sub$true.froh)
    mean.called <- mean(sub$call.froh)
    mean.abs.diff <- mean(sub$abs.diff)
    sd.abs.diff <- sd(sub$abs.diff)
    mean.diff <- mean(sub$froh.diff)
    sd.diff <- sd(sub$froh.diff)
    save <- c(mean.true, mean.called, mean.abs.diff, sd.abs.diff, mean.diff, sd.diff, unlist(sub[1,c(5:14)]))
    OUT1 <- rbind(OUT1, save)
  } else{
    next
  }
  ct <- ct+1
}
colnames(OUT1)[1:6] <- c('mean.true.froh','mean.called.froh','mean.abs.diff','sd.abs.diff','mean.diff','sd.diff')
pop.lev.res <- as.data.frame(OUT1)
for(c in 1:6){
  pop.lev.res[,c] <- as.numeric(pop.lev.res[,c])
}
pop.lev.res$covg <- factor(pop.lev.res$covg, levels = c(5, 10, 15, 30, 50))
pop.lev.res$phws <- factor(pop.lev.res$phws, levels = c(50, 100, 1000))

# write.csv(OUT1[,c(1,2,6:ncol(OUT1))], 'population_froh_results_SA.csv', row.names = FALSE)
# write.csv(froh.stats[,c(1:3,6:ncol(froh.stats))], 'individual_froh_results_SA.csv', row.names = FALSE)