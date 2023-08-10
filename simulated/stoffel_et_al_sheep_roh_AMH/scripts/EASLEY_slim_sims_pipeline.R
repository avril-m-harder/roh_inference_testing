library(data.table)
#library(tidyverse)
source("/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/scripts/make_slim.R")
source("/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/scripts/combine_mut_roh.R")
source("/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/scripts/correct_vcf.R")
library(furrr)
library(future)
library(glue)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

#### FOR TESTING #####
# demo_mod <- 'TEST'
# pop_size1 <- 100
# pop_size2 <- 100
# pop_size3 <- 100
# pop_size4 <- 50
# time1 <- 1000
# time2 <- 2000
# time3 <- 3000
# time4 <- 4000
######################

# run slim simulation
args_in <- commandArgs(trailingOnly=TRUE)
print(args_in)
if (length(args_in) != 9) stop("Provide demo model name, pop_size1:4, and time1:4")
demo_mod <- args_in[1]
pop_size1 <- as.numeric(args_in[2])
pop_size2 <- as.numeric(args_in[3])
pop_size3 <- as.numeric(args_in[4])
pop_size4 <- as.numeric(args_in[5])
time1 <- as.numeric(args_in[6])
time2 <- as.numeric(args_in[7])
time3 <- as.numeric(args_in[8])
time4 <- as.numeric(args_in[9])

# directory structure
out_dir <- paste0('/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/',
                  demo_mod)
# to_create <- paste0('/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/',
#                     demo_mod, "/", c("muts", "out", "roh", "slim_code", "trees", "vcfs"))
# walk(to_create, dir.create, recursive = TRUE, showWarnings = TRUE)

# dotdotdot is input for making the slim file, see make_slim.R
slim_roh <- function(seeds, pop_size = pop_size1, ...) {
   
   run_name <- paste0("tasdev_", seeds)
   
   # create slim_simulation
   make_slim(seed = seeds, out_dir = out_dir, ...)
   
   # run slim
   #system(paste0("slim -time -Memhist -seed ", seed, " slim_sim/sims/slim_code/tasdev_", seed, ".slim"))
   system(paste0("slim -time -seed ", seeds," ", out_dir, "/slim_code/tasdev_", demo_mod, "_", seeds, ".slim"))
   
   # recapitation and overlay of neutral mutations
   # check that folders are there

   system(paste("python /scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/scripts/slim2_overlay_mut.py", run_name, demo_mod, pop_size1)) 
   
   # correct vcf 
   correct_vcf(paste0("/scratch/avrilh/roh_param_project/roh_inference_testing/simulated/stoffel_et_al_sheep_roh_AMH/output/", demo_mod, "/vcfs/tasdev_", seeds, ".vcf"))

}


## parameters
mut1_dom_coeff <- 0.05   ## previously c(0, 0.05, 0.2)
mut1_gam_mean <- -0.03   ## previously c(-0.01, -0.03, -0.05)
mut1_gam_shape <- 0.2
genome_size <- 30e6               ## in Stoffel script, 1e8 (100 Mb)
mut_rate_del <- (1.39e-8 * 0.7)   ## in Stoffel script, 70% deleterious = 1e-8 * 0.7 = 7e-9
recomb_rate <- 9.15e-9            ## uniform

params <- data.table::transpose(as.data.frame(c(pop_size1, pop_size2, pop_size3, pop_size4, mut1_dom_coeff, mut1_gam_mean,
                                                mut1_gam_shape, genome_size, mut_rate_del, recomb_rate, time1, time2, time3, time4)))
colnames(params) <- c("pop_size1", "pop_size2", "pop_size3", "pop_size4", "mut1_dom_coeff", "mut1_gam_mean", "mut1_gam_shape",
                      "genome_size", "mut_rate_del", "recomb_rate", "time1", "time2", "time3", "time4")

# params <- c(pop_size1, pop_size2, pop_size3, pop_size4, mut1_dom_coeff, mut1_gam_mean, mut1_gam_shape,
#             genome_size, mut_rate_del, recomb_rate, time1, time2, time3, time4) %>%
#    setNames(c("pop_size1", "pop_size2", "pop_size3", "pop_size4", "mut1_dom_coeff", "mut1_gam_mean", "mut1_gam_shape",
#               "genome_size", "mut_rate_del", "recomb_rate", "time1", "time2", "time3", "time4")) %>%
#    as_tibble()

## only running 1 replicate for now (and possibly forever)
num_sim_per_parset <- 1
set.seed(123)
seeds <- sample(1:1e5, num_sim_per_parset * nrow(params))
# replicate each parameter set num_sim_per_parset times
params_sim <- params[rep(1:nrow(params), each = num_sim_per_parset), ] %>%
   mutate(seed = seeds)

# plan(multiprocess, workers = 32) ## AMH: commented this out because tested locally 
# make all roh and trees files and recapitate
future_pmap(params_sim, slim_roh, pop_size1)
