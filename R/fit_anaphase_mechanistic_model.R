  library(bayesplot)
  library(tidybayes)
  library(ggplot2)
  library(patchwork)

fit_anaphase_mechanistic_model <- function(jobset_str, t_ana_df, K=160,
                                        identifier = 'anaphase_multiple_v211',
                                        run_analysis = TRUE,
                                        use_parallel = FALSE,
					fits_folder_str = 'fits',
					plot_opt = 1,
					num_iter = 1000,
					load_existing = FALSE,
					stan_file = here::here('src/stan_files/anaphase_reversals_hierarchical.stan'),
					use_3D_model = TRUE
){
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
  identifier = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  Data <- process_jobset(jobset_str,K = K,max_missing = 0.25,plot_opt = plot_opt)
  if (is.infinite(K)) { K = max(Data$Frame) } #find final frame and use this if defaulting to all the data
  dt = Data$Time %>% unique() %>% diff() %>% min()
  pairIDs <- unique(Data$SisterPairID)
  y = prepare_for_stan_format(Data)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  } #replace missing values with zero for stan compatibility
  y_missing <- purrr::map(y_missing,as.integer)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  run_mcmc_for_single_trajectory <- function(trajectory_index,Data,dt,K,y,y_missing,pairIDs,t_ana_df,
					     run_analysis=TRUE,load_existing=FALSE){
  nTracks <- 1
  if (load_existing && file.exists(file.path(fits_folder_str,paste('anaphase_mechanistic_sister_',pairIDs[trajectory_index],identifier,'.rds',sep='')))){
    estimate = readRDS(file.path(fits_folder_str,paste('anaphase_mechanistic_sister_',pairIDs[trajectory_index],identifier,'.rds',sep='')))
  } else {
  if (run_analysis){
    cos_phi <- get_cos_phi(Data,pairIDs[trajectory_index])
    cos_phi[is.na(cos_phi)] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
    stan_input = list(dt=dt, T=K,
                      nStates = 5,
                      y = y[[trajectory_index]],
                      sigma0 = c(0,0.5,0.5,0,0),
		      T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0),
		      T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1),
		      t_ana_input = t_ana_df[["t_ana"]][trajectory_index],
		      cos_phi = cos_phi
    )
    initF <- function() list(tau = rep(450,nTracks)+100*rnorm(nTracks),
                             v_plus = rep(0.03,nTracks), v_minus = rep(-0.03,nTracks),
                             alpha = rep(0.01,nTracks), kappa = rep(0.05,nTracks),
                             p_coh=0.9+0.02*rnorm(nTracks), p_icoh=0.9+0.02*rnorm(nTracks),
                             v_ana = rep(0.04,nTracks), t_ana = t_ana_df[["t_ana"]][trajectory_index]
    )
    if (use_3D_model){
      cat("Warning: overwriting chosen stan model with anaphase_3d.stan instead due to option use_3D=TRUE \n")
      stan_file = here::here('src/stan_files/anaphase_3d.stan')
    }
    estimate <- rstan::stan(file=stan_file,
                     data=stan_input,
                     seed = 42,
                     chains = 4,
                     warmup = num_iter,
                     iter = 2*num_iter,
		     pars = c("eta","f","auxStates","P","p_ana","aux"),
		     include=FALSE, #avoid saving the params listed above
                     control=list(adapt_delta=0.95,
                                  max_treedepth=12))
    saveRDS(estimate, file = file.path(fits_folder_str,paste('anaphase_mechanistic_sister_',pairIDs[trajectory_index],identifier,'.rds',sep='')))
  while ((!check_max_Rhat_less_than_tol(estimate)$converged) & num_iter < 100){
    #check if converged using default tolerance value
    #if not, will need to run again
    num_iter <- num_iter*2
    cat(paste("Run for trajectory ", trajectory_index, " did not converge. Trying this trajectory again with ", num_iter, "iterations instead."))
        estimate <- rstan::stan(file=stan_file,
                     data=stan_input,
                     seed = 42,
                     chains = 4,
                     warmup = num_iter,
                     iter = 2*num_iter,
                     pars = c("eta","f","auxStates","P","p_ana","aux"),
                     include=FALSE, #avoid saving the params listed above
                     control=list(adapt_delta=0.95,
                                  max_treedepth=12))
  }
  } else {
    #just load from previous save
    estimate = readRDS(file.path(fits_folder_str,paste('anaphase_mechanistic_sister_',pairIDs[trajectory_index],identifier,'.rds',sep='')))
  }
  }
    if (plot_opt){
      plot_anaphase_mechanistic_results_single_trajectory(jobset_str,estimate,Data,dt,pairIDs[trajectory_index],identifier)
    }
  return(estimate)
  }
  if (use_parallel){
    library(furrr)
    plan(multiprocess,workers=parallel::detectCores()) #on nero need to specify how many workers here
    my_map <- furrr::future_map
  } else {
    my_map <- purrr::map
  }
anaphase_estimates <- my_map(seq_along(pairIDs),
                             function(x) try(run_mcmc_for_single_trajectory(x,Data,dt,K,y,y_missing,pairIDs,t_ana_df,run_analysis=run_analysis)))
  return(anaphase_estimates)
}

plot_anaphase_mechanistic_results <- function(jobset_str,estimate,Data,dt,pairIDs,
                                           show_mcmc_diagnostics=FALSE){
  draws <- estimate %>%
    spread_draws(t_ana[track]) %>%
    mutate(SisterPairID = pairIDs[track])

  p1 <- draws %>%
    ggplot(aes(y = SisterPairID, x = t_ana)) +
    geom_eyeh() +
    theme_bw() +
    labs(x="Time of anaphase transition (s)",
         y = "Sister Pair ID")

  D_sim_draws <- estimate %>% spread_draws(D_sim[track,timept]) %>%
    summarize(lb = quantile(D_sim, probs = 0.025),
              median = quantile(D_sim, probs = 0.5),
              ub = quantile(D_sim, probs = 0.975)) %>%
    mutate(Time = timept*dt,
           SisterPairID = pairIDs[track])

  t_ana_inferred <- estimate %>% spread_draws(t_ana[track]) %>%
    summarize(lb = quantile(t_ana, probs = 0.025),
              median = quantile(t_ana, probs = 0.5),
              ub = quantile(t_ana, probs = 0.975)) %>%
    mutate(SisterPairID = pairIDs[track])

  p2 <- Data %>% group_by(SisterPairID,Time) %>%
    summarise(intersister_dist = abs(diff(Position_1))) %>% left_join(D_sim_draws) %>%
    ggplot(aes(x=Time,y=intersister_dist)) +
    geom_line() +
    facet_wrap(.~SisterPairID) +
    geom_ribbon(aes(x=timept*dt,ymax=ub,ymin=lb),alpha=0.3) +
    geom_vline(data=t_ana_inferred,aes(xintercept=lb),color="red",linetype="dotted") +
    geom_vline(data=t_ana_inferred,aes(xintercept=median),color="red") +
    geom_vline(data=t_ana_inferred,aes(xintercept=ub),color="red",linetype="dotted") +
    labs(y="Intersister distance (um)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90))
  print(p1+p2)
#  ggsave(stringr::str_replace(jobset_str,".csv",
#                              paste("_anaphase_results_",identifier,".eps",sep="")),device=cairo_ps)
  if (show_mcmc_diagnostics){
    #plots for mcmc diagnostics
    pars_to_plot = paste(c("tau","alpha","beta","a","t_ana"),"[",rep(1:nTracks,1,each=5),"]",sep="")
    color_scheme_set("purple")
    posterior <- as.array(estimate)
    mcmc_hist(posterior,pars=pars_to_plot)
    mcmc_dens_overlay(posterior,pars=pars_to_plot)
    mcmc_trace(posterior,pars=pars_to_plot)
  }
  return (p1+p2)
}

fit_anaphase_mechanistic_model_single_trajectory <- function(jobset_str, t_ana_input, K=160,
                                           identifier = 'anaphase_single_v211',
                                           run_analysis = TRUE,
                                           trajectory_index = 1,
					   fits_folder_str='fits',
					   plot_opt = 1,
					   num_iter = 1000,
					   stan_file = here::here('src/stan_files/anaphase_reversals_hierarchical.stan'),
					   use_3D_model = TRUE,
                                           dt = NULL
){
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
  identifier = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  Data <- process_jobset(jobset_str,K = K,max_missing = 0.25,plot_opt = plot_opt) %>%
    filter(!is.na(SisterID)) #remove unpaired KT trajectories
  if (is.null(dt)) { dt=Data$Time %>% unique() %>% diff() %>% min() }
  pairIDs <- unique(Data$SisterPairID)
  nTracks = 1
  y = prepare_for_stan_format(Data)
  y_missing = map(y, is.na) %>% map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  } #replace missing values with zero for stan compatibility
  y_missing <- purrr::map(y_missing,as.integer)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  if (run_analysis){
    cos_phi <- get_cos_phi(Data,pairIDs[trajectory_index])
    cos_phi[is.na(cos_phi)] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
    stan_input = list(dt=dt, T=K,
                      nStates = 5,
                      y = y[[trajectory_index]],
                      y_missing = y_missing[[trajectory_index]],
                      sigma0 = c(0,0.5,0.5,0,0),
                      t_ana_input=t_ana_input,
                      T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0),
                      T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1),
		      cos_phi = cos_phi
    )
    pos_part <- function(q) max(q,0)
    initF <- function() list(tau = rep(450,nTracks)+100*rnorm(nTracks),
                             v_plus = sapply(rep(0.03,nTracks) + 0.02*rnorm(nTracks),pos_part), 
			     v_minus = sapply(rep(-0.03,nTracks) + 0.02*rnorm(nTracks),pos_part),
                             alpha = sapply(rep(0.01,nTracks) + 0.01*rnorm(nTracks),pos_part), 
			     kappa = sapply(rep(0.03,nTracks) + 0.02*rnorm(nTracks),pos_part),
                             p_coh = sapply(0.9+0.1*rnorm(nTracks),pos_part), 
			     p_icoh = sapply(0.9+0.1*rnorm(nTracks),pos_part),
                             v_ana = sapply(rep(0.03,nTracks)+0.02*rnorm(nTracks),pos_part), 
			     t_ana = t_ana_input+10*dt*rnorm(nTracks),
			     gamma = 0.05
    )
    if (use_3D_model){
      cat("Warning: overwriting chosen stan model with anaphase_3d.stan instead due to option use_3D=TRUE \n")
      stan_file = here::here('src/stan_files/anaphase_3d.stan')
    }
    m <- stan_model(stan_file)
    estimate <- sampling(m,
                              data=stan_input,
                              seed = 42,
                              chains = 4,
	                      warmup = num_iter,
         	              iter = 2*num_iter,
                              pars = c("eta","f","auxStates","P","p_ana","aux"),
                              include=FALSE, #avoid saving the params listed above
                              control=list(adapt_delta=0.95,
                                           max_treedepth=12))
    saveRDS(estimate, file = file.path(fits_folder_str,paste('anaphase_multiple_estimates',identifier,'.rds',sep='')))
  } else {
    #just load from previous save
    estimate = readRDS(file.path(fits_folder_str,paste('anaphase_multiple_estimates',identifier,'.rds',sep='')))
  }
  if (plot_opt){
    plot_anaphase_mechanistic_results_single_trajectory(jobset_str,estimate,Data,dt,pairIDs[trajectory_index],identifier)
  }
  return(estimate)
}

plot_anaphase_mechanistic_results_dissolving_spring_single_trajectory <- function(
    jobset_str,estimate,Data,dt,pairID,identifier,
    show_mcmc_diagnostics=FALSE){

  draws <- estimate %>%
    spread_draws(t_ana,t_spring) %>%
    mutate(SisterPairID = pairID)

  p1 <- draws %>%
    tidyr::gather(event,time,-.chain,-.iteration,-.draw,-SisterPairID) %>%
    ggplot(aes(y = event, x = time)) +
    geom_eyeh(aes(color=event)) +
    theme_bw() +
    labs(x="Time of anaphase transition (s)",
         y = "Event")

  D_sim_draws <- estimate %>% spread_draws(D_sim[timept]) %>%
    summarize(lb = quantile(D_sim, probs = 0.025),
              median = quantile(D_sim, probs = 0.5),
              ub = quantile(D_sim, probs = 0.975)) %>%
    mutate(Time = timept*dt,
           SisterPairID = pairID)

  t_ana_inferred <- estimate %>% spread_draws(t_ana,t_spring) %>%
    summarize(lb = quantile(t_ana, probs = 0.025),
              median = quantile(t_ana, probs = 0.5),
              ub = quantile(t_ana, probs = 0.975),
              lb_spring = quantile(t_spring, probs = 0.025),
              median_spring = quantile(t_spring, probs = 0.5),
              ub_spring = quantile(t_spring, probs = 0.975)) %>%
    mutate(SisterPairID = pairID)

  p2 <- Data %>% group_by(SisterPairID,Time) %>%
    filter(SisterPairID == pairID) %>%
    summarise(intersister_dist = abs(diff(Position_1))) %>% left_join(D_sim_draws) %>%
    ggplot(aes(x=Time,y=intersister_dist)) +
    geom_line() +
    facet_wrap(.~SisterPairID) +
    geom_ribbon(aes(x=timept*dt,ymax=ub,ymin=lb),alpha=0.3) +
    geom_vline(data=t_ana_inferred,aes(xintercept=lb),color="red",linetype="dotted") +
    geom_vline(data=t_ana_inferred,aes(xintercept=median),color="red") +
    geom_vline(data=t_ana_inferred,aes(xintercept=ub),color="red",linetype="dotted") +
    geom_vline(data=t_ana_inferred,aes(xintercept=lb_spring),color="cyan",linetype="dotted") +
    geom_vline(data=t_ana_inferred,aes(xintercept=median_spring),color="cyan") +
    geom_vline(data=t_ana_inferred,aes(xintercept=ub_spring),color="cyan",linetype="dotted") +
    labs(y="Intersister distance (um)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90))
  print(p1+p2)
#  ggsave(stringr::str_replace(jobset_str,".csv",
#                              paste("_anaphase_results_",identifier,".eps",sep="")),device=cairo_ps)
  if (show_mcmc_diagnostics){
    #plots for mcmc diagnostics
    #TODO: this would require fixing
    pars_to_plot = paste(c("tau","alpha","beta","a","t_ana"),"[",rep(1:nTracks,1,each=5),"]",sep="")
    color_scheme_set("purple")
    posterior <- as.array(changept_estimate)
    mcmc_hist(posterior,pars=pars_to_plot)
    mcmc_dens_overlay(posterior,pars=pars_to_plot)
    mcmc_trace(posterior,pars=pars_to_plot)
  }
  return (p1+p2)
}

plot_anaphase_mechanistic_results_single_trajectory <- function(
  jobset_str,estimate,Data,dt,pairID,identifier,
  show_mcmc_diagnostics=FALSE){

  draws <- estimate %>%
    spread_draws(t_ana) %>%
    mutate(SisterPairID = pairID)

  p1 <- draws %>%
    ggplot(aes(y = SisterPairID, x = t_ana)) +
    geom_eyeh() +
    theme_bw() +
    labs(x="Time of anaphase transition (s)",
         y = "Sister Pair ID")

  D_sim_draws <- estimate %>% spread_draws(D_sim[timept]) %>%
    summarize(lb = quantile(D_sim, probs = 0.025),
              median = quantile(D_sim, probs = 0.5),
              ub = quantile(D_sim, probs = 0.975)) %>%
    mutate(Time = timept*dt,
           SisterPairID = pairID,
	   Frame=timept)

  t_ana_inferred <- estimate %>% spread_draws(t_ana) %>%
    summarize(lb = quantile(t_ana, probs = 0.025),
              median = quantile(t_ana, probs = 0.5),
              ub = quantile(t_ana, probs = 0.975)) %>%
    mutate(SisterPairID = pairID)

  p2 <- Data %>% group_by(SisterPairID,Frame) %>%
    filter(SisterPairID == pairID) %>%
    summarise(intersister_dist = abs(diff(Position_1))) %>%
    left_join(D_sim_draws,by=c("SisterPairID","Frame")) %>%
    ggplot(aes(x=(Frame-1)*dt,y=intersister_dist)) +
    geom_line() +
    facet_wrap(.~SisterPairID) +
    geom_ribbon(aes(x=(Frame-1)*dt,ymax=ub,ymin=lb),alpha=0.3) +
    geom_vline(data=t_ana_inferred,aes(xintercept=lb),color="red",linetype="dotted") +
    geom_vline(data=t_ana_inferred,aes(xintercept=median),color="red") +
    geom_vline(data=t_ana_inferred,aes(xintercept=ub),color="red",linetype="dotted") +
    labs(y="Intersister distance (um)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90))
  print(p1+p2)
  #  ggsave(stringr::str_replace(jobset_str,".csv",
  #                              paste("_anaphase_results_",identifier,".eps",sep="")),device=cairo_ps)
  if (show_mcmc_diagnostics){
    #plots for mcmc diagnostics
    pars_to_plot = paste(c("tau","alpha","beta","a","t_ana"),"[",rep(1:nTracks,1,each=5),"]",sep="")
    color_scheme_set("purple")
    posterior <- as.array(changept_estimate)
    mcmc_hist(posterior,pars=pars_to_plot)
    mcmc_dens_overlay(posterior,pars=pars_to_plot)
    mcmc_trace(posterior,pars=pars_to_plot)
  }
  return (p1+p2)
}


run_anaphase_mechanistic_analysis_for_folder <- function(path_to_folder,
                                           K_list = NA,
                                           identifier = "JHv001",
                                           run_analysis = TRUE,
					   use_parallel = FALSE,
					   fits_folder_str = 'fits')
{
  jobset_str_list = list.files(path = path_to_folder,pattern="\\.csv$",full.names=TRUE)
  if (is.na(K_list)){ #set default
    K_list <- rep(Inf,length(jobset_str_list))
  } else if (length(K_list)==1) { #if integer input rather than vector
    K_list <- rep(K_list,length(jobset_str_list))
  }
  out <- purrr::map2(jobset_str_list, K_list, function(x,y) fit_anaphase_mechanistic_model(x, K=y,
                                        identifier = identifier,
                                        run_analysis = run_analysis,
                                        use_parallel = use_parallel,
					fits_folder_str = fits_folder_str))
  return(out)
}

plot_anaphase_mechanistic_batch_results <- function(out,jobset_str,K=Inf,dt=4.5,plot_opt=0){

  Data <- process_jobset(jobset_str,K = K,max_missing = 0.25,plot_opt = plot_opt)
  pairIDs <- unique(Data$SisterPairID)
  draws <- out %>% map(function(x) spread_draws(x,t_ana)) %>%
    bind_rows(.id = "track") %>%
    mutate(track = as.integer(track)) %>%
    mutate(SisterPairID = pairIDs[track])
  p1 <- draws %>%
    ggplot(aes(y = SisterPairID, x = t_ana)) +
    geom_eyeh() +
    theme_bw() +
    labs(x="Time of anaphase transition (s)",
         y = "Sister Pair ID")
  print(p1)

  D_sim_draws <- out %>% map(function(x) spread_draws(x,D_sim[timept]) %>%
    summarize(lb = quantile(D_sim, probs = 0.025),
              median = quantile(D_sim, probs = 0.5),
              ub = quantile(D_sim, probs = 0.975))) %>%
    bind_rows(.id="track") %>%
    mutate(track = as.integer(track)) %>%
    mutate(SisterPairID = pairIDs[track],
           Time = timept*dt)

t_ana_inferred <- draws %>% group_by(SisterPairID) %>%
    summarize(lb = quantile(t_ana, probs = 0.025),
              median = quantile(t_ana, probs = 0.5),
              ub = quantile(t_ana, probs = 0.975))

 p2 <- Data %>% group_by(SisterPairID,Time) %>%
    summarize(intersister_dist = abs(diff(Position_1))) %>% left_join(D_sim_draws) %>%
    ggplot(aes(x=Time,y=intersister_dist)) +
    geom_line() +
    facet_wrap(.~SisterPairID) +
    geom_ribbon(aes(x=timept*dt,ymax=ub,ymin=lb),alpha=0.3) +
    geom_vline(data=t_ana_inferred,aes(xintercept=lb),color="red",linetype="dotted") +
    geom_vline(data=t_ana_inferred,aes(xintercept=median),color="red") +
    geom_vline(data=t_ana_inferred,aes(xintercept=ub),color="red",linetype="dotted") +
    labs(y="Intersister distance (um)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90))
  print(p1+p2)
}

