setwd('/Users/Avril/Documents/roh_param_project/roh_inference_testing/simulated/data/plink_output/')
fns <- list.files(path = 'plink_round2/', pattern = 'hom')

## select which round to process for sensitivity analysis
rnd <- 'round2'

##### Read in PLINK results and summarize #####
## Loop over demographic scenarios and write 1 file for each
demos <- unique(do.call(rbind, strsplit(fns, split = '_'))[,1])

for(d in demos){
  print(d)
  called.roh.id.counter <- 1
  
  OUT <- matrix(c('id','start','end','n.snps','covg','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','called.roh.id'), nrow = 1)
  write.table(OUT, paste0(rnd,'/',d,'_PLINK_all_coordinates.txt'),
              quote = FALSE, row.names = FALSE, sep='\t', col.names = FALSE)
  
  d.fns <- fns[grep(d, fns)]
  for(f in d.fns){
    vars <- unlist(strsplit(f, split = '_'))
    vars[22] <- gsub('.hom','',vars[22])
    phwh <- vars[8]
    phwm <- vars[10]
    phws <- vars[12]
    phzd <- vars[14]
    phzg <- vars[16]
    phwt <- vars[18]
    phzs <- vars[20]
    phzk <- vars[22]
    covg <- vars[4]
    dat <- read.table(paste0('plink_',rnd,'/',f), header=TRUE)
    s.id <- called.roh.id.counter
    e.id <- called.roh.id.counter + nrow(dat) - 1 
    if(nrow(dat) > 0){
      dat <- dat[,c(1,7,8,10)]
      dat$covg <- covg
      dat$phwh <- phwh
      dat$phwm <- phwm
      dat$phws <- phws
      dat$phzd <- phzd
      dat$phzg <- phzg
      dat$phwt <- phwt
      dat$phzs <- phzs
      dat$phzk <- phzk
      dat$called.roh.id <- c(s.id:e.id)

      write.table(dat, paste0(rnd,'/',d,'_PLINK_all_coordinates.txt'),
                  append = TRUE, quote = FALSE, row.names = FALSE, sep='\t', col.names = FALSE)
      
      called.roh.id.counter <- called.roh.id.counter + nrow(dat)
    }
  }
}
