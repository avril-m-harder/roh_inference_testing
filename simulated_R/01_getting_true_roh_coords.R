setwd('/Users/Avril/Desktop/roh_param_project/simulated_plotting_summarizing/bcftoolsroh_output/')
library(ape)
library(scales)
`%notin%` <- Negate(`%in%`)

##### Read in SLiM results #####
## set chromosome length for SLiM run
c.len <- 30e6

## read in ancestral fasta:
## homozygous, so read in twice.
## bases not included in VCF assumed to still be homozygous.
anc.1 <- read.FASTA('/Users/Avril/Desktop/roh_param_project/SLiM/22_01_06_testing_slim/ancestral.fasta')
anc.1 <- as.character(anc.1)

## read in VCF
vcf <- read.table('/Users/Avril/Desktop/roh_param_project/SLiM/22_01_06_testing_slim/final_pop.vcf', header=TRUE)
head(vcf[,c(1:10)])

## convert reference alleles in VCF to nuc at that position in ancestral fasta
bases <- c('A','T','C','G')
for(i in 1:nrow(vcf)){
  ref <- toupper(anc.1$`1`[vcf$POS[i]])
  vcf[i, 'REF'] <- ref
  vcf[i, 'ALT'] <- sample(bases[bases %notin% ref], 1)
}
all(vcf$ALT != vcf$REF) ## verify that all REF and ALT alleles differ for each locus (should be TRUE)

##### Check ROH / heterozygosity statistics #####
## check ROH length distributions based on variants in VCF
samp.cols <- grep('i', colnames(vcf)) ## get columns that contain sample data

OUT <- NULL    ## where to save het sites for each individual
# OUT1 <- NULL   ## where to save true ROH coordinates for each individual
for(col in samp.cols){
  print(colnames(vcf)[col])
  samp <- gsub('i','',colnames(vcf)[col])
  temp <- vcf[,c(2,col)] ## subset to sample genotypes
  temp[temp[,2] == '0|1', 2] <- '1'
  temp[temp[,2] == '1|0', 2] <- '1'
  temp[temp[,2] == '1|1', 2] <- '0'
  temp[temp[,2] == '0|0', 2] <- '0'
  het.sites <- temp[temp[,2] == 1,]
  het.sites <- het.sites[order(het.sites$POS),]
  het.sites[,2] <- samp
  colnames(het.sites) <- c('het.pos','samp')
  OUT <- rbind(OUT, het.sites)
  
  # for(h in 1:nrow(het.sites)){
  #   if(h < nrow(het.sites) & het.sites[h, 1] == (het.sites[h+1, 1] + 1)){
  #     print(paste0(colnames(vcf)[col]," - CONSECUTIVE HET SITES"))
  #   } else{
  #     if(h == 1){
  #       s <- 1
  #       e <- het.sites[h, 1] - 1
  #       save <- c(samp, s, e, (e-s+1))
  #       OUT1 <- rbind(OUT1, save)
  #     } 
  #     if(h == nrow(het.sites)){
  #       s <- het.sites[h-1, 1] + 1
  #       e <- het.sites[h, 1] - 1
  #       save <- c(samp, s, e, (e-s+1))
  #       OUT1 <- rbind(OUT1, save)
  #       s <- het.sites[h, 1] + 1
  #       e <- c.len
  #       save <- c(samp, s, e, (e-s+1))
  #       OUT1 <- rbind(OUT1, save)
  #     } 
  #     if(h != 1 & h != nrow(het.sites)){
  #       s <- het.sites[h-1, 1] + 1
  #       e <- het.sites[h, 1] - 1
  #       save <- c(samp, s, e, (e-s+1))
  #       OUT1 <- rbind(OUT1, save)
  #     }
  #   }
  # }
  if(col == min(samp.cols)){
    write.table(OUT, '../true_indiv_data/true_het_sites.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
    # write.table(OUT1, '../true_indiv_data/true_roh_coords.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
    OUT <- NULL
    # OUT1 <- NULL
  } else{
    write.table(OUT, '../true_indiv_data/true_het_sites.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
    # write.table(OUT1, '../true_indiv_data/true_roh_coords.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
    OUT <- NULL
    # OUT1 <- NULL
  }
}



