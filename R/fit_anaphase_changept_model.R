fit_anaphase_changept_model <- function(jobset_str, K=Inf,
                                        identifier = 'changept_multiple_v111',
                                        run_analysis = TRUE,
					fits_folder_str=here::here('fits'),
					plot_opt = 1,
					max_missing = 0.25,
					dt = NULL,
					num_iter = 1000
){
  job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
  identifier = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  Data <- process_jobset(jobset_str,K = K,max_missing = max_missing,plot_opt = plot_opt) %>%
        filter(!is.na(SisterID))
  stopifnot(nrow(Data)>0) # ensure there is some data to fit the model to for this job
  if (is.infinite(K)) { K = max(Data$Frame) } #find final frame and use this if defaulting to all the data
  if (is.null(dt)){ dt = Data$Time %>% unique() %>% diff() %>% min() }

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
  nTracks = length(pairIDs)
  stopifnot(nTracks>10) #ensure at least 10 pairs to fit to
  y = prepare_for_stan_format(Data)
  y_missing = map(y, is.na) %>% map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  } #replace missing values with zero for stan compatibility
  y_missing <- map(y_missing,as.integer)
  if (run_analysis){
    stan_input = list(dt=dt, T=K,
                      nTracks = nTracks,
                      y = y,
                      y_missing = y_missing
    )
    changept_estimate <- stan(file=here::here('src/stan_files/anaphase_changepoint_multiple_tracks.stan'),
                              data=stan_input,
                              seed = 42,
                              chains = 4,
	                      warmup = num_iter,
         	              iter = 2*num_iter,
                              control=list(adapt_delta=0.95,
                                           max_treedepth=12))
    saveRDS(changept_estimate, file = file.path(fits_folder_str,paste('anaphase_changepoint_estimates',identifier,'.rds',sep='')))
  } else {
    #just load from previous save
    changept_estimate = readRDS(file.path(fits_folder_str,paste('anaphase_changepoint_estimates',identifier,'.rds',sep='')))
  }
  if (plot_opt){
    plot_anaphase_changept_results(jobset_str,changept_estimate,Data,dt,pairIDs)
  }
  return(changept_estimate)
}

plot_anaphase_changept_results <- function(jobset_str,changept_estimate,Data,dt,pairIDs,
                                           show_mcmc_diagnostics=FALSE){
  library(bayesplot)
  library(tidybayes)
  library(ggplot2)
  library(patchwork)
  draws <- changept_estimate %>%
    spread_draws(t_ana[track]) %>%
    mutate(SisterPairID = pairIDs[track])

  p1 <- draws %>%
    ggplot(aes(y = SisterPairID, x = t_ana)) +
    geom_eyeh() +
    theme_bw() +
    labs(x="Time of anaphase transition (s)",
         y = "Sister Pair ID")

  D_sim_draws <- changept_estimate %>% spread_draws(D_sim[track,timept]) %>%
    summarize(lb = quantile(D_sim, probs = 0.025),
              median = quantile(D_sim, probs = 0.5),
              ub = quantile(D_sim, probs = 0.975)) %>%
    mutate(Time = timept*dt,
           SisterPairID = pairIDs[track])

  t_ana_inferred <- changept_estimate %>% spread_draws(t_ana[track]) %>%
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
  ggsave(stringr::str_replace(jobset_str,".csv",
                              paste("_anaphase_results_",identifier,".eps",sep="")))
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

