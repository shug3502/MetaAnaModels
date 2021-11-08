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
source(here::here('R/fit_anaphase_changept_model.R'))
source(here::here('R/fit_anaphase_reversals_model.R'))
source(here::here('R/make_tracking_figure.R'))
source(here::here('R/make_population_params_and_heterogeneity_figure.R'))
source(here::here('R/make_anaphase_times_and_speed_figure.R'))
source(here::here('R/ana_spatial_analysis.R'))
source(here::here('R/make_force_profiles_and_maturation_figure.R'))
source(here::here('R/generate_figures_based_on_states_and_switching.R'))
source(here::here('R/switching_analysis_helper_fns.R'))
source(here::here('R/hawkes_process_analysis.R'))
source(here::here('R/make_coordination_of_anaphase_onset_figure.R'))
source(here::here('R/make_local_coordination_agreement_figure.R'))
source(here::here('R/get_directional_switch_events.R'))
source(here::here('R/extract_hidden_states.R'))
##############################
#options etc
dt=2.05
fits_folder_str <- "fits"
identifier = "jonathanharrison_MetaAnaModels_2s_v411" #args[1]
path_to_folder = "data" #args[2]

jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",
                              full.names=TRUE,recursive=TRUE)
K_list <- rep(Inf,length(jobset_str_list))
stopifnot(length(jobset_str_list)>0)

draws <- summarise_batch_anaphase_reversals(identifier,jobset_str_list,dt=dt,
                                     fits_folder_str = fits_folder_str,
                                     max_missing = 0.25,tol = 1.05,
                                     run_analysis = FALSE,
                                     use_parallel = FALSE)
#############################

#for Figure 1: tracking
make_tracking_figure(jobset_str_list)

#############################
#Figures 2 to 4 generated separately or to be added later
#############################
#Figure 5: Distribution of population parameters and heterogeneity
make_population_params_and_heterogeneity_figure(jobset_str_list,draws)

#############################
#Figure 6: Force profiles and maturation
make_force_profile_and_maturation_figure(jobset_str_list,draws,identifier,
                                         num_iter=200,dt=dt,min_num_sisters=10,nframes_early=300,
                                         fits_folder_str=fits_folder_str)

#############################
#Figure 7: anaphase times and speeds
make_anaphase_times_and_speed_figure(jobset_str_list,draws,min_num_sisters=10,nframes_early=300)

# #############################
# #Figure 8: Directional switching
# #TODO: sort this out by looping over jobs or changin how it works
# generate_figures_based_on_states_and_switching(estimate,sigma_sim,jobset_str,identifier,dt=dt)
# 
# #############################
# #Figure 9: Coordination of switches via Hawkes process
# #TODO: similarly
# hawkes_process_analysis(sigma_sim,jobset_str,dt)
# 
# #############################
# #Figure 10: local coordination and agreement
# #TODO: similarly
# make_local_coordination_agreement_figure(jobset_str,estimate,dt=dt,
#                                                      nStates=6,niter=200)
# 
# #############################
 #Figure 11: coordination of anaphase onset
make_coordination_of_anaphase_onset_figure(identifier,jobset_str_list,dt=dt)

for (i in seq_along(jobset_str_list)){
  #extract hidden states
#  jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame//kittracking001-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture10_flowdec_deconvolved.ome.csv")
  job_id = stringr::str_split(jobset_str_list[i],"kittracking")[[1]][2]
  edited_job_id = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  path_to_est <- file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',
                                                 edited_job_id,'.rds',sep=''))
  if (file.exists(path_to_est)){
    estimate <- readRDS(file=path_to_est)
    sigma_sim <- extract_hidden_states(estimate)
#     success = generate_figures_based_on_states_and_switching(estimate,sigma_sim,jobset_str_list[i],
# 							     paste0(identifier,"cell",i),dt)
    #############################
    #Figure 8: Directional switching
    #TODO: sort this out by looping over jobs or changin how it works
    generate_figures_based_on_states_and_switching(estimate,sigma_sim,jobset_str_list[i],
                                                   paste0(identifier,"cell",i),dt=dt)
    
    #############################
    #Figure 9: Coordination of switches via Hawkes process
    #TODO: similarly
    hawkes_process_analysis(sigma_sim,jobset_str_list[i],dt)
    
    #############################
    #Figure 10: local coordination and agreement
    #TODO: similarly
    make_local_coordination_agreement_figure(jobset_str_list[i],estimate,dt=dt,
                                             nStates=6,niter=200)
    
  } else {
    cat('Files from sampling do not exist. Skipping this cell.\n')
    success=0
  }
  if (success){cat("SUCCESS!\n")} else {cat("FAILED :(\n")}
}
