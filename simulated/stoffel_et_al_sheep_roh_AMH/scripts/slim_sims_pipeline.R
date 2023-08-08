library(data.table)
#library(tidyverse)
source("make_slim.R")
source("combine_mut_roh.R")
source("correct_vcf.R")
library(furrr)
library(future)
library(glue)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

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

#### FOR TESTING #####
# demo_mod <- 'TEST'
# pop_size1 <- 1000
# pop_size2 <- 750
# pop_size3 <- 500
# pop_size4 <- 250
# time1 <- 1000
# time2 <- 2000
# time3 <- 3000
# time4 <- 4000
######################


# directory structure
out_dir <- paste0('/Users/Avril/Documents/roh_param_project/manuscript_resubmission/revised_SLiM/stoffel_et_al_sheep_roh-master_AMH/output/',
                  demo_mod)
to_create <- paste0('/Users/Avril/Documents/roh_param_project/manuscript_resubmission/revised_SLiM/stoffel_et_al_sheep_roh-master_AMH/output/',
                    demo_mod, "/", c("muts", "out", "roh", "slim_code", "trees", "vcfs"))
walk(to_create, dir.create, recursive = TRUE, showWarnings = TRUE)

# dotdotdot is input for making the slim file, see make_slim.R
slim_roh <- function(seeds, pop_size = pop_size1, ...) {
   
   # run_name <- paste0("tasdev_", demo_mod, "_", seeds)
   run_name <- paste0("tasdev_", seeds)
   
   # create slim_simulation
   make_slim(seed = seeds, out_dir = out_dir, ...)
   
   # run slim
   #system(paste0("slim -time -Memhist -seed ", seed, " slim_sim/sims/slim_code/tasdev_", seed, ".slim"))
   system(paste0("slim -time -seed ", seeds," ", out_dir, "/slim_code/tasdev_", demo_mod, "_", seeds, ".slim"))
   
   # recapitation and overlay of neutral mutations
   # check that folders are there

   system(paste("python3 slim2_overlay_mut.py", run_name, demo_mod, pop_size1)) 
   
   # correct vcf (AMH: commented out on Aug. 8, 2023 - correcting REF/ALT alleles appears to be unnecessary)
   # correct_vcf(paste0("../output/", demo_mod, "/vcfs/tasdev_", seeds, ".vcf"))

   ##### AMH: still going to make nucleotide-specific VCF files and simulate seq reads from there #####   
   # call ROH
   # use vcf output to call ROH
   # currently, pyslim sometimes outputs weird REF and ALT positions in the vcf file
   # these might be the mutations from SLiM
   # When one of these is the first line in the VCF gt table, then REF is the empty
   # string and ALT is a large number. plink doesn't like this and will find no ROH
   # This means that some of the simulations can't be processed further.
   # system(paste0("plink --vcf ", out_path, "/vcfs/", run_name, ".vcf ",  # /usr/local/bin/plink
   #               "--tasdev --out ", out_path, "/roh/", run_name, " ",
   #               "--homozyg --homozyg-window-snp 25 --homozyg-snp 25 --homozyg-kb 390 ",
   #               "--homozyg-gap 250 --homozyg-density 100 --homozyg-window-missing 0 ",
   #               "--homozyg-het 0 ",
   #               "--homozyg-window-het 0"))
   
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



##### Original code - not using for now because 1 combo & 1 rep per demo scenario #####
# # combinations
# mut1_dom_coeff <- 0.05   ## previously c(0, 0.05, 0.2)
# mut1_gam_mean <- -0.03   ## previously c(-0.01, -0.03, -0.05)
# mut1_gam_shape <- 0.2
# genome_size <- 30e6               ## in Stoffel script, 1e8 (100 Mb)
# mut_rate_del <- (1.39e-8 * 0.7)   ## in Stoffel script, 70% deleterious = 1e-8 * 0.7 = 7e-9
# recomb_rate <- 9.15e-9            ## uniform
# 
# params <- c(pop_size1, pop_size2, pop_size3, pop_size4, mut1_dom_coeff, mut1_gam_mean, mut1_gam_shape,
#                       genome_size, mut_rate_del, recomb_rate) %>%
#    setNames(c("pop_size1", "pop_size2", "pop_size3", "pop_size4", "mut1_dom_coeff", "mut1_gam_mean", "mut1_gam_shape",
#               "genome_size", "mut_rate_del", "recomb_rate")) %>%
#    as_tibble() %>%
#    mutate(time1 = time1,        
#           time2 = time2,
#           time3 = time3,
#           time4 = time4)
# 
# ## only running 1 replicate for now (and possibly forever)
# num_sim_per_parset <- 1
# set.seed(123)
# seeds <- sample(1:1e5, num_sim_per_parset * nrow(params))
# # replicate each parameter set num_sim_per_parset times
# params_sim <- params[rep(1:nrow(params), each =num_sim_per_parset), ] %>%
#    mutate(seed = seeds)
# 
# # plan(multiprocess, workers = 32) ## AMH: commented this out because tested locally 
# # make all roh and trees files and recapitate
# future_pmap(params_sim, slim_roh, pop_size1)
# 
# # make all roh and trees files and recapitate
# # future_map(seeds, slim_roh,  genome_size = 1e8, pop_size1 = 1000, pop_size2 = 200, 
# #                   time1 = 1000, time2 = 1100, mut_rate_del = 1e-9, recomb_rate = 1e-8, # 1.27e-8
# #                   # mutationtype 1
# #                   mut1_dom_coeff = 0.05, mut1_gam_mean = -0.03, mut1_gam_shape = 0.2)
# #                   # mutationtype2
# #                   # mut2_dom_coeff = 0.01, mut2_gam_mean = -0.1, 
# #                   # mut2_gam_shape = 1.5, mut2_rel_freq = 0.1)
# #                   #mut3_gam_mean = 0.03, 
# #                   #mut3_gam_shape = 2, mut3_rel_freq = 0.01, 
# #                   #mut3_dom_coeff = 0.8
# #                   #)
# 
# # make a safe combine function
# combine_safe <- safely(combine_mut_roh)
# # combine mutations and roh data and calculate length classes
# out <- map(paste0("tasdev_", seeds), combine_safe, out_path, 
#            roh_cutoff_small = 1560, roh_cutoff_long = 12500)
# # roh_cutoff_small = 1560, roh_cutoff_long = 6250)
# 
# # extract only non-error runs and combine
# mut_df <- out %>% 
#    purrr::transpose() %>% 
#    .[[1]] %>% 
#    bind_rows() %>% 
#    mutate(id = paste0(ind_id, seed),
#           seed = as.numeric(seed)) %>% 
#    dplyr::arrange(id) %>% 
#    left_join(params_sim, by = "seed")
# # complete(nesting(id, roh_class)) %>% 
# # arrange(id)
# 
# #saveRDS(mut_df, file = "slim_sim/sims/out/sims_weakly_03.RData")
# 
# write_delim(mut_df, paste0(out_path, "/out/par_combs_popsize1_", pop_size1,
#                            "_popsize2_", pop_size2, "_popsize3_", pop_size3, ".txt"))
# 
# 
# # mut_p <- mut_df %>% 
# #    group_by(seed, roh_class) %>% 
# #    summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE))
# # 
# # source("../tasdev_ID/theme_simple.R")
# # library(gghalves)
# # 
# # p1 <- ggplot(mut_p, aes(roh_class, s_sum_per_MB, fill = roh_class)) +
# #    geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
# #                    transformation_params = list(height = 0, width = 1.3, seed = 1)) +
# #    geom_half_boxplot(side = "r", outlier.color = NA,
# #                      width = 0.6, lwd = 0.3, color = "black",
# #                      alpha = 0.8) +
# #    #scale_y_continuous(limits = c(-0.02, 0)) +
# #    scale_fill_viridis_d() +
# #    ylab("selection coefficient per cM") + 
# #  #  ggtitle("Sim: weakly del\nROH cutoff 4900KB, \nmut1. dist: -0.03, 2, dom coeff 0.1 \nmut2. dist: -0.2, 3, dom coeff 0.01") +
# #    # geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
# #    theme_simple(grid_lines = FALSE, axis_lines = TRUE)
# # p1
# # ggsave("figs/del_roh_dom005_sel003.jpg", plot = p1, width = 6, height = 4)
# # 
# # mut_df
# # 
# # library(lme4)
# # 
# # mod <- lmer(s_sum_per_MB ~ roh_class + (1|seed), data = mut_df)
# # summary(mod)
# # tidy(mod, conf.int = TRUE)
# # library(broom.mixed)
# # 
# # tidy(mod, conf.int = TRUE)
