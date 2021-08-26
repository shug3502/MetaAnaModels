args = commandArgs(trailingOnly=TRUE)
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
source(here::here('R/hawkes_process_analysis.R'))
source(here::here('R/get_directional_switch_events.R'))
source(here::here('R/extract_hidden_states.R'))
source(here::here('R/make_tracking_figure.R'))
##############################
#options etc
dt=2.05
fits_folder_str <- "fits"
identifier = args[1]
path_to_folder = args[2]
#############################

#for Figure 1
make_tracking_figure()

jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",full.names=TRUE)
K_list <- rep(Inf,length(jobset_str_list))
stopifnot(length(jobset_str_list)>0)

for (i in seq_along(jobset_str_list)){
  #extract hidden states
#  jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame//kittracking001-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture10_flowdec_deconvolved.ome.csv")
  job_id = stringr::str_split(jobset_str_list[i],"kittracking")[[1]][2]
  edited_job_id = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  path_to_est <- file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',edited_job_id,'.rds',sep=''))
#  path_to_est <- here::here("../Constandia/estimate.rds")
  if (file.exists(path_to_est)){
    estimate <- readRDS(file=path_to_est)
    sigma_sim <- extract_hidden_states(estimate)
#    sigma_sim <- readRDS(here::here("../Constandia/sigma_sim.rds"))
    success = generate_figures_based_on_states_and_switching(estimate,sigma_sim,jobset_str_list[i],
							     paste0(identifier,"cell",i),dt)
  } else {
    cat('Files from sampling do not exist. Skipping this cell.\n')
    success=0
  }
  if (success){cat("SUCCESS!\n")} else {cat("FAILED :(\n")}

hawkes_fit <- hawkes_process_analysis(sigma_sim,jobset_str_list[i],dt)
}
