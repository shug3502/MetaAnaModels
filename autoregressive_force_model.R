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
source('R/helper_fns.R')
##############################
#options etc
run_analysis <- TRUE
num_iter <- 1000
fits_folder_str <- "fits"
dt=2.05
identifier = "force_dcomposition"
path_to_folder = "data"
set.seed(123)
# jobset_str <- sample(jobset_str_list,1)
jobset_str <- "data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking006-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture5_flowdec_deconvolved.ome.csv"
Data_2s <- process_jobset(jobset_str,max_missing=0.25,K=100,
                          plot_opt=FALSE) %>%
  dplyr::mutate(filename=jobset_str)

plot_output <- function(fit,pid){
  s <- summary(fit,pars="mt_force")$summary
  tt <- rownames(s) %>% stringr::str_extract(pattern="\\d{1,3},") %>% stringr::str_extract("\\d{1,3}") %>% as.integer()
  id <- rep(c(1,2),max(tt))
  
  df <- tibble(f = s[,"50%"],
               frame=tt,
               SisterID=id)
  p1 <- ggplot(df,aes((frame-1)*dt,f,color=factor(SisterID))) + 
    geom_line() + 
    theme_minimal() +
    scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(1,4)]) +
    labs(x="Time (s)",
         y="Force (um/s)",
         color="SisterID")
  
  p2 <- df %>% tidyr::spread(key=SisterID,value=f) %>% 
    ggplot(aes(`1`,`2`)) + 
    geom_point() +
    theme_minimal() +
    labs(x="Force on sister 1 (um/s)",
         y="Force on sister 2 (um/s)")
  
  p1 | p2
  ggsave(here::here(paste0("plots/force_decomposition_",identifier,"_pair",pid,".eps")),
         width=6,height=4)
}

fit_model <- function (pid){
 dat <- Data_2s %>% 
    filter(SisterPairID==pid)
  K <- max(dat$Frame)
  y = prepare_for_stan_format(dat)
  y_missing = map(y, is.na) %>% map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  } #replace missing values with zero for stan compatibility
  y_missing <- purrr::map(y_missing,as.integer)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  if (run_analysis){
    cos_phi <- get_cos_phi(dat,pid)
    cos_phi[is.na(cos_phi)] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
    stan_input = list(dt=dt, T=K,
                      y = y[[1]],
                      y_missing = y_missing[[1]],
                      T0 = start_end$T0,
                      T1 = start_end$T1,
                      cos_phi = cos_phi
    )
    pos_part <- function(q) abs(q)
    nTracks <- 1
    initF <- function() list(tau = rep(450,nTracks)+100*rnorm(nTracks),
                             alpha = sapply(rep(0.01,nTracks) + 0.01*rnorm(nTracks),pos_part), 
                             kappa = sapply(rep(0.03,nTracks) + 0.02*rnorm(nTracks),pos_part)
                             )
    stan_file <- here::here("src/stan_files/metaphase_decommpose_MT_force.stan")
    m <- stan_model(stan_file)
    estimate <- sampling(m,
                         data=stan_input,
                         seed = 42,
                         chains = 4,
                         warmup = num_iter,
                         iter = 2*num_iter,
                         init=initF,
                         control = list(adapt_delta=0.95)
                         )
    saveRDS(estimate, file = file.path(fits_folder_str,paste('metaphase_',identifier,'force_decomposition',pid,'.rds',sep='')))
  } else {
    #just load from previous save
    estimate = readRDS(file.path(fits_folder_str,paste('metaphase_',identifier,'force_decomposition',pid,'.rds',sep='')))
  }
  plot_output(estimate,pid)
  return(estimate)
}

pair_ids <- unique(Data_2s$SisterPairID)
out_ana <- purrr::map(pair_ids,fit_model)

fit <- fit_model(10)

get_num_negative_frames <- function(fit){
  s <- summary(fit,pars="mt_force")$summary
  tt <- rownames(s) %>% stringr::str_extract(pattern="\\d{1,3},") %>% stringr::str_extract("\\d{1,3}") %>% as.integer()
  id <- rep(c(1,2),max(tt))
  
  df <- tibble(f = s[,"50%"],
               frame=tt,
               SisterID=id)
  
  df %>% group_by(frame) %>%
    summarise(all_negative = all(f<0)) %>%
    ungroup() %>%
    summarise(num_neg_frames = sum(all_negative)) %>% #max(sequence(rle(all_negative)$lengths)*(all_negative))) %>%
    pull(num_neg_frames)
}

approximate_discrete_states <- function(fit){
s <- summary(fit,pars="mt_force")$summary
tt <- rownames(s) %>% stringr::str_extract(pattern="\\d{1,3},") %>% stringr::str_extract("\\d{1,3}") %>% as.integer()
id <- rep(c(1,2),max(tt))
df <- tibble(f = s[,"50%"],
             frame=tt,
             SisterID=id)
df %>% arrange(SisterID) %>% group_by(frame) %>%
  summarise(sigma = case_when((first(f)>0) & (last(f)>0) ~ "++",
                              (first(f)>0) & (last(f)<0) ~ "+-",
                              (first(f)<0) & (last(f)>0) ~ "-+",
                              (first(f)<0) & (last(f)<0) ~ "--")) %>%
  ggplot(aes(frame,sigma)) + 
  geom_point() + theme_minimal()
}

n <- out_ana %>% purrr::map_dbl(get_num_negative_frames)
tibble(SisterPairID=pair_ids,n=n) %>%
  inner_join(u) %>%
  ggplot(aes(any_usable,n)) + 
  # geom_violin(draw_quantiles=0.5) +
  # geom_jitter() + 
  geom_boxplot() +
  theme_minimal()
#######################

check_rhat <- function(object){
  max(summary(object,pars=c("alpha","kappa","L","tau",
                            "mt_force","phi","epsilon"))$summary[,"Rhat"])
}
r <- purrr::map_dbl(out_ana,check_rhat)
d <- purrr::map_dbl(out_ana,rstan::get_num_divergent)

tibble(SisterPairID=pair_ids,n=n) %>%
       inner_join(u) %>% inner_join(tibble(SisterPairID = pair_ids,
                                    rhat = r,
                                    divergences = d) %>%
                                      mutate(converged = case_when(is.na(rhat) ~ "missing",
                                                            rhat > 1.1 ~ "unconverged",
                                                            divergences > 2 ~ "divergences",
                                                            rhat < 1.1 ~ "converged",
                                                            TRUE ~ "unknown"))) %>%
  mutate(might_help = (!(converged=="converged") & !any_usable)) %>% pull(might_help)


