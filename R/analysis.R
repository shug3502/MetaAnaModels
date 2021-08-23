library(rstan) 
library(dplyr) 
library(purrr)
library(zoo) 
library(here) 
library(ggplot2) 
library(patchwork) 
library(bayesplot)
library(tidybayes)
library(stringr)
library(readr)
library(tidyr)

rstan::rstan_options(auto_write = TRUE) #tries to avoid recompiling stan code
options(mc.cores = parallel::detectCores()) #uses as many cores as you have
source(here::here('R/helper_fns.R'))
source(here::here('R/switching_analysis_helper_fns.R'))
source(here::here('R/extract_hidden_states.R'))
##############################
#options etc
dt=2.05
identifier <- "JH_test_v000"
fits_folder_str <- "fits"
#############################

  #extract hidden states
  jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame//kittracking001-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture10_flowdec_deconvolved.ome.csv")
  job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
  edited_job_id = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  path_to_est <- file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',edited_job_id,'.rds',sep=''))
  if (file.exists(path_to_est)){
    estimate <- readRDS(file=path_to_est)
    sigma_sim <- extract_hidden_states(estimate)
    success = generate_figures_based_on_states_and_switching(estimate,sigma_sim,jobset_str,
							     paste0(identifier,"cell",1),dt)
  } else {
    cat('Files from sampling do not exist. Skipping this cell.\n')
    success=0
  }
  if (success){cat("SUCCESS!\n")} else {cat("FAILED :(\n")}



