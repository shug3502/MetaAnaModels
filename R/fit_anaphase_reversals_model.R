fit_anaphase_reversals_model <- function(jobset_str, t_ana_df, K=Inf,
                                        identifier = 'anaphase_reversals_v111',
                                        run_analysis = TRUE,
                                        use_parallel = FALSE,
					fits_folder_str = here::here('fits'),
					plot_opt = 1,
					num_iter = 1000,
					max_missing=0.25,
					dt = NULL,
				        stan_file = here::here('src/stan_files/anaphase_reversals_hierarchical.stan')
){
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
  identifier = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  Data <- process_jobset(jobset_str,K = K,max_missing = max_missing,plot_opt = plot_opt) %>%
	filter(!is.na(SisterID))
  if (is.infinite(K)) { K = max(Data$Frame) } #find final frame and use this if defaulting to all the data
  if (is.null(dt)) { dt=Data$Time %>% unique() %>% diff() %>% min() }

  #model for anaphase can have problems if data only contains metaphase so ensure chromatids for each pair separate
  #use a cut off on the 1D kk distance at the last tracked frame to do this
  has_separated_df <- Data %>%
    group_by(SisterPairID,Frame) %>%
    arrange(SisterID,SisterPairID,Frame) %>%
    summarise(kk_dist=-diff(Position_1)) %>%                
    group_by(SisterPairID) %>% 
    summarise(kk_dist_last = last(kk_dist[!is.na(kk_dist)]), 
              separated = kk_dist_last>2.0) #cut off of 2um  
  Data <- Data %>% inner_join(has_separated_df %>% ungroup %>% filter(separated),
                              by="SisterPairID")

  pairIDs <- unique(Data$SisterPairID)
  nTracks <- length(pairIDs)
  stopifnot(nTracks>10) #ensure at least 10 pairs to fit to
  y = prepare_for_stan_format(Data)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  } #replace missing values with zero for stan compatibility
  y_missing <- purrr::map(y_missing,as.integer)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(SisterPairID=pairIDs[as.integer(traj_index)]) #use this info to only fit to the existing tracks

  if (run_analysis){
    cos_phi <- purrr::map(pairIDs,function(x) get_cos_phi(Data,x))
    for (ii in seq_along(cos_phi)){
      cos_phi[[ii]][is.na(cos_phi[[ii]])] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
    }
    stan_input = list(dt=dt, T=K,
                      y = y,
                      sigma0 = c(0,0.5,0.5,0,0),
		      T0 = start_end %>% pull(T0),
		      T1 = start_end %>% pull(T1),
		      t_ana_input = t_ana_df[["t_ana"]],
		      cos_phi = cos_phi
    )
    initF <- function() list(tau = rep(450,nTracks)+100*rnorm(nTracks),
                             v_plus = rep(0.03,nTracks), v_minus = rep(-0.03,nTracks),
                             alpha = rep(0.01,nTracks), kappa = rep(0.05,nTracks),
                             p_coh=0.9+0.02*rnorm(nTracks), p_icoh=0.9+0.02*rnorm(nTracks),
                             v_ana = rep(0.04,nTracks), t_ana = t_ana_df[["t_ana"]] + 3*dt*rnorm(nTracks)
    )
    #stan_file = here::here('src/stan_files/anaphase_reversals_hierarchical.stan')
    m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
    estimate <- rstan::sampling(m,
                     data=stan_input,
                     seed = 42,
                     chains = 4,
                     warmup = num_iter,
                     iter = 2*num_iter,
		     pars = c("eta","f","auxStates","P","p_ana","aux","state_probs","conditional_state_probs","P_col"),
		     include=FALSE, #avoid saving the params listed above
                     control=list(adapt_delta=0.95,
                                  max_treedepth=15))
    saveRDS(estimate, file = file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',identifier,'.rds',sep='')))
  while ((!check_max_Rhat_less_than_tol(estimate)$converged) & num_iter < 100){
    #check if converged using default tolerance value 
    #if not, will need to run again
    num_iter <- num_iter*2
    cat(paste("Run did not converge. Trying again with ", num_iter, "iterations instead."))
    estimate <- rstan::sampling(m,
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
    estimate = readRDS(file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',identifier,'.rds',sep='')))
  }
  if (plot_opt){
    plot_anaphase_mechanistic_results(jobset_str,estimate,Data,dt,pairIDs)
  }
  return(estimate)
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

