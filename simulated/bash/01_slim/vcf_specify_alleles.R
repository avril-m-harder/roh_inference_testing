library(ape)
library(readr)

`%notin%` <- Negate(`%in%`)

setwd('/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/data/01_slim/output')

## read in ancestral fasta:
## bases not included in VCF assumed to still be homozygous.
anc.1 <- read.FASTA('ancestral.fasta')
anc.1 <- as.character(anc.1)

## read in VCF
fns <- list.files(pattern = "final*")

for(f in fns){
  demo <- strsplit(f, split = '_')[[1]][2]
  new.fn <- paste0(demo,'_n50_spec_allele.vcf')
  vcf <- read.table(f, header=FALSE, skip = 6)
  colnames(vcf) <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',paste0('i',c(1:50)))
  
  ## remove triallelic loci (also including a fourth allele filter, just in case that shows up in really large pops)
  TO.RM <- NULL
  for(r in 1:nrow(vcf)){
    if(length(grep("2", vcf[r, grep("i", colnames(vcf))])) > 0){
      TO.RM <- c(TO.RM, r)
    }
    if(length(grep("3", vcf[r, grep("i", colnames(vcf))])) > 0){
      TO.RM <- c(TO.RM, r)
    }
  }
  if(length(TO.RM) > 0){
	  print(paste0(f,'\n\n'))
	  print(vcf[TO.RM,])
	  vcf <- vcf[-TO.RM,]
  }
  
  ## convert reference alleles in VCF to nuc at that position in ancestral fasta
  bases <- c('A','T','C','G')
  for(i in 1:nrow(vcf)){
    ref <- toupper(anc.1$`1`[vcf$POS[i]])
    ## account for possibility of undefined alleles in reference FASTA file
    if(ref == 'N'){
      ref <- sample(bases, 1)
    }
    vcf[i, 'REF'] <- ref
    vcf[i, 'ALT'] <- sample(bases[bases %notin% ref], 1)
  }
  print(paste0(demo,' - allele ID check: ',all(vcf$ALT != vcf$REF))) ## verify that all REF and ALT alleles differ for each locus (should be TRUE)
  
  meta <- readr::read_delim(f, delim = "\t", n_max = 5, 
                            col_names = FALSE, show_col_types = FALSE)[[1]]
  readr::write_lines(meta, file = new.fn)
  # append genotypes to make file complete
  data.table::fwrite(vcf, file = new.fn, append = TRUE, 
                     col.names = TRUE, sep = "\t")
}
