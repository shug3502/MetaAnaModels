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
#library(furrr)
rstan::rstan_options(auto_write = TRUE) #tries to avoid recompiling stan code
options(mc.cores = parallel::detectCores()) #uses as many cores as you have
#plan(multisession, workers = 16)
source('R/helper_fns.R')
source('R/switching_analysis_helper_fns.R')
##############################
#options etc
run_analysis <- TRUE
is_zero_force_model <- TRUE
num_iter <- 2000
fits_folder_str <- "fits"
subsample_factor <- 1
dt = 2.05*subsample_factor
path_to_folder = args[1]
identifier = if_else(is_zero_force_model,paste0("understanding_model_misspecification_zero_force_",args[2]),
                     paste0("understanding_model_misspecification_4states_",args[2])) #args[1]
jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",
                              full.names=TRUE,recursive=TRUE)
# stopifnot(length(jobset_str_list)>0)
set.seed(123)
# jobset_str <- sample(jobset_str_list,1)
#jobset_str <- "data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking006-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture5_flowdec_deconvolved.ome.csv"
for (jobset_str in jobset_str_list){
jobid = jobset_str %>% stringr::str_extract(pattern="kittracking.+_flowdec_deconvolved") #kittracking part of string
print(jobid)
Data_2s <- process_jobset(jobset_str,max_missing=0.25,K=Inf,
                          plot_opt=FALSE,interpolation=FALSE) %>%
#  dplyr::filter((Frame %% subsample_factor) == 1) %>%
  dplyr::mutate(filename=jobset_str)

fit_model <- function (pid){
  cat("fitting moodel for pair ", pid, ".... \n")
  dat <- Data_2s %>% 
    filter(SisterPairID==pid)
  K <- dat$Frame %>% unique %>% length
  y = prepare_for_stan_format(dat)
  y_missing = map(y, is.na) %>% map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  } #replace missing values with zero for stan compatibility
  y_missing <- purrr::map(y_missing,as.integer)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  if (run_analysis){
if (!file.exists(file.path(fits_folder_str,paste('anaphase_',jobid,identifier,'poor_oscillating_track',pid,'.rds',sep='')))){
    cos_phi <- get_cos_phi(dat,pid)
    cos_phi[is.na(cos_phi)] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
    nStates = if_else(is_zero_force_model,5,4)
    stan_input = list(dt=dt, T=K,
                      nStates = nStates,
                      y = y[[1]],
                      y_missing = y_missing[[1]],
                      sigma0 = rep(1,nStates)/nStates,
                      T0 = start_end$T0,
                      T1 = start_end$T1,
                      cos_phi = cos_phi
    )
    pos_part <- function(q) abs(q)
    nTracks <- 1
    changept_file <- here::here("src/stan_files/anaphase_changepoint.stan")
    changept_model <- stan_model(changept_file)
    changept_fit <- sampling(changept_model,
                         data=stan_input,
                         seed = 42,
                         chains = 4,
                         warmup = num_iter,
                         iter = 2*num_iter)
    t_ana_est <- rstan::extract(changept_fit,pars="t_ana")$t_ana %>% median()
    stan_input$t_ana_input <- t_ana_est    
    initF <- function() list(tau = rep(450,nTracks)+100*rnorm(nTracks),
                             v_plus = sapply(rep(0.03,nTracks) + 0.02*rnorm(nTracks),pos_part), 
                             v_minus = -sapply(rep(0.03,nTracks) + 0.02*rnorm(nTracks),pos_part),
                             alpha = sapply(rep(0.01,nTracks) + 0.01*rnorm(nTracks),pos_part), 
                             kappa = sapply(rep(0.03,nTracks) + 0.02*rnorm(nTracks),pos_part),
                             v_ana = sapply(rep(0.03,nTracks) + 0.02*rnorm(nTracks),pos_part),
                             t_ana = sapply(rep(t_ana_est,nTracks) + 5*dt*rnorm(nTracks),pos_part)
    )
    stan_file <- here::here("src/stan_files/anaphase_missing.stan")
    m <- stan_model(stan_file)
    estimate <- sampling(m,
                         data=stan_input,
                         seed = 42,
                         chains = 4,
                         warmup = num_iter,
                         iter = 2*num_iter,
                         pars = c("eta","f","auxStates","P","p_ana","aux"),
                         init=initF,
                         include=FALSE, #avoid saving the params listed above
                         control=list(adapt_delta=0.95,
                                      max_treedepth=12))
    if (identical(as.array(estimate),numeric(0))){
      cat("MCMC output is empty\n")
      cat(paste0("t_A estimated via changepoint model as ",t_ana_est, " compared to total movie length ", K*dt,"\n"))
    }
    saveRDS(estimate, file = file.path(fits_folder_str,paste('anaphase_',jobid,identifier,'poor_oscillating_track',pid,'.rds',sep='')))
} else {
  cat("file ", file.path(fits_folder_str,paste('anaphase_',jobid,identifier,'poor_oscillating_track',pid,'.rds',sep='')), " already exists\n")
  estimate <- NULL
}
  } else {
    #just load from previous save
    estimate = readRDS(file.path(fits_folder_str,paste('anaphase_',jobid,identifier,'poor_oscillating_track',pid,'.rds',sep='')))
  }
  return(estimate)
}

pair_ids <- unique(Data_2s$SisterPairID)
#out_ana <- furrr::future_map(pair_ids,fit_model)
out_ana <- purrr::map(pair_ids,fit_model)
}
# 


