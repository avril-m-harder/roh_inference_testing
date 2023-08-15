library(vcfR)

setwd('/Users/Avril/Desktop/')

vcf <- read.vcfR('decline_sample_cvg_50x_filtered_SNPs.vcf')
gts <- extract.gt(vcf)

for(c in 1:ncol(gts)){
  gts[,c] <- gsub('0/0','0', gts[,c], fixed = TRUE)
  gts[,c] <- gsub('1/1','0', gts[,c], fixed = TRUE)
  gts[,c] <- gsub('1|1','0', gts[,c], fixed = TRUE)
  
  gts[,c] <- gsub('0/1','1', gts[,c], fixed = TRUE)
  gts[,c] <- gsub('0|1','1', gts[,c], fixed = TRUE)
  print(  table(gts[,c]))
}

OUT <- NULL
for(c in 1:ncol(gts)){
  save <- c(colnames(gts)[c], length(gts[which(gts[,c] == 1), c]), length(gts[which(gts[,c] == 1), c])/nrow(gts))
  OUT <- rbind(OUT, save)
}

hist(as.numeric(OUT[,3]))
OUT[which(as.numeric(OUT[,3]) < 0.05),]

## removing samples i24 and i38