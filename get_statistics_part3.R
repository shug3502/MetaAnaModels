library(naniar)
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
library(furrr)
plan(multisession, workers = 16)
rstan::rstan_options(auto_write = TRUE) #tries to avoid recompiling stan code
options(mc.cores = parallel::detectCores()) #uses as many cores as you have
source('R/helper_fns.R')
source('R/switching_analysis_helper_fns.R')

path_to_folder = "data"
is_zero_force_model <- TRUE
run_id <- "v332"
identifier = if_else(is_zero_force_model,paste0("understanding_model_misspecification_zero_force_","jonathanharrison_",run_id),
                     paste0("understanding_model_misspecification_4states_","jonathanharrison_",run_id))
subsample_factor = 1
fits_folder_str <- "fits"

jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",
                              full.names=TRUE,recursive=TRUE)


################
check_rhat <- function(object){
  max(summary(object,pars="theta")$summary[,"Rhat"])
}
check_any_samples <- function(object){
  is_nonempty <- !(identical(as.array(object),numeric(0)))
  return(is_nonempty)
}

get_diagnostics_df <- function(path_to_fit){
  kittracking = path_to_fit %>% stringr::str_extract(pattern="kittracking.+_flowdec_deconvolved")
  pair = path_to_fit %>% stringr::str_extract(pattern="_track\\d{1,3}") %>%
    stringr::str_extract(pattern="\\d{1,3}") %>% as.integer()
  fit <- readRDS(path_to_fit)
  is_nonempty = check_any_samples(fit)
  if (is_nonempty){
    r = check_rhat(fit)
    d = rstan::get_num_divergent(fit)
    nchains = length(dimnames(as.array(fit))$chains)
  } else {
    r = NA
    d = NA
    nchains = 0
  }
  return(tibble(kittracking=kittracking,
                SisterPairID = pair,
                is_nonempty = is_nonempty,
                rhat = r,
                nchains = nchains,
                divergences = d))
}
#########################

get_median_params <- function(path_to_fit){
  fit <- readRDS(path_to_fit)
  is_nonempty <- check_any_samples(fit)
  if (is_nonempty){
    out <- tibble(param=c("tau","alpha","kappa","v_minus","v_plus","p_icoh","p_coh","L","v_ana","t_ana"),
           theta=rstan::summary(fit,pars="theta")$summary[,"50%"])
  } else {
    out <- tibble(param=c("tau","alpha","kappa","v_minus","v_plus","p_icoh","p_coh","L","v_ana","t_ana"),
         theta=rep(NA,10))
  }
  return(out)
}

get_results <- function(jobset_str){
  Data_2s <- process_jobset(jobset_str,max_missing=0.25,K=Inf,
                            plot_opt=FALSE,interpolation=FALSE) %>%
 #   dplyr::filter((Frame %% subsample_factor) == 1) %>%
    dplyr::mutate(filename=jobset_str)
  pair_ids <- unique(Data_2s$SisterPairID)
  jobid = jobset_str %>% stringr::str_extract(pattern="kittracking.+_flowdec_deconvolved") #kittracking part of string

  get_median_params_track <- function(pid){
    path_to_fit <- file.path(fits_folder_str,paste('anaphase_',jobid,identifier,'poor_oscillating_track',pid,'.rds',sep=''))
    out_df <- get_median_params(path_to_fit)
    return(out_df)
  }
if (length(pair_ids)==0){
  res <- tibble(SisterPairID = rep(NA,10), param=c("tau","alpha","kappa","v_minus","v_plus","p_icoh","p_coh","L","v_ana","t_ana"),
         theta=rep(NA,10))
} else {
  res <- purrr::map_df(pair_ids,get_median_params_track,.id="SisterPairID") %>% mutate(SisterPairID = pair_ids[as.integer(SisterPairID)])
}
  return(res)
}
#########################
draws <- purrr::map_df(jobset_str_list,get_results,.id="cell") %>% 
  mutate(filename=jobset_str_list[as.integer(cell)]) %>%
  mutate(p_reversal=NA,p_revtoana=NA) %>%
  as_tibble()

pro_fun <- function(jobset_str,start_from=0,max_missing=0.25,K=Inf) read.csv(jobset_str,header=TRUE) %>% 
  as_tibble() %>%
  dplyr::ungroup() %>%
  extract_long_tracks(K,start_from,max_missing) %>%
  reorder_sisters() %>%
  exclude_sisters_that_cross()
Data <- purrr::map(jobset_str_list,pro_fun) %>% 
  bind_rows(.id="cell") %>%
  mutate(filename = jobset_str_list[as.integer(cell)])

ana_times_df <- draws %>% 
  filter(param=="t_ana") %>%
  group_by(filename) %>%
  summarise(t_ana = median(theta,na.rm=T))

intersister_dist_df <- Data %>% 
  left_join(ana_times_df,by="filename") %>%
  mutate(metaphase = (Time>=0) & (Time < t_ana)) %>%
  filter(metaphase) %>%
  group_by(filename,SisterPairID,Frame) %>%
  summarise(intersister_dist = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2),
            radius = mean(sqrt(Position_2^2+Position_3^2),na.rm=T),
            centroid_x = mean(Position_1,na.rm=T)) %>%
  group_by(filename,SisterPairID) %>%
  summarise(av_intersister_dist = median(intersister_dist,na.rm=TRUE),
            av_radius = median(radius,na.rm=T),
            av_x = mean(centroid_x,na.rm=T),
            DAP = mean(abs(centroid_x-av_x),na.rm=T))

missing_df <- Data %>%
  group_by(filename,SisterPairID,SisterID) %>%
  summarise(num_missing = sum(is.na(Position_1)),
            proportion_missing=sum(is.na(Position_1))/length(Position_1)) %>%
  group_by(filename,SisterPairID) %>%
  summarise(proportion_missing = max(proportion_missing,na.rm=T),
	    num_missing = max(num_missing,na.rm=T))

#check for convergence for each Sister Pair individually
path_to_fits <- list.files("fits",pattern=paste0("jonathanharrison_",run_id),full.names=T)
diagnostics_df <- furrr::future_map(path_to_fits,get_diagnostics_df) %>%
  bind_rows() %>%
  mutate(model="anaphase_missing")

df <- left_join(draws,intersister_dist_df,by = c("SisterPairID", "filename")) %>% 
  left_join(missing_df, by=c("filename","SisterPairID")) %>%
  mutate(kittracking = stringr::str_extract(filename,
	 pattern="kittracking.+_flowdec_deconvolved")) %>%
  left_join(diagnostics_df, by = c("SisterPairID","kittracking"))


saveRDS(df,paste0("fits/draws_median_anaphase_",run_id,".rds"))

df2 <- df %>% group_by(filename,SisterPairID) %>%
    summarise(converged = all(rhat<1.1),
	      unconverged = any(rhat>1.1),
	      failed = sum(is.na(rhat))) %>% 
    group_by(filename) %>%
    summarise(nconverged = sum(converged,na.rm=T),
              nunconverged=sum(unconverged,na.rm=T),
	      nfailed=sum(failed,na.rm=T)) 
p1 <- ggplot(df2,aes(nconverged,nunconverged)) + 
    geom_point() +
    theme_bw()
p2 <- ggplot(df2,aes(nconverged,nfailed)) +
    geom_point() +
    theme_bw()
p1 | p2
  ggsave(here::here(paste0("plots/conv_vs_unconv_",run_id,".eps")),height=7,width=7)

#files1 = draws$filename %>% unique()
#files2 = intersister_dist_df$filename %>% unique()
#files3 = missing_df$filename %>% unique()
#files4 = diagnostics_df$kittracking %>% unique()
#setdiff(files1,files2)
#setdiff(files1,files,3)
#stringr::str_extract(files1,pattern="kittracking.+_flowdec_deconvolved") %>% setdiff(.,files4)

#for (pp in path_to_fits){
#print(pp)
#a <- get_diagnostics_df(pp)
#}
