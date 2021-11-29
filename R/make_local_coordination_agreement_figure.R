make_local_coordination_agreement_figure <- function(jobset_str_list,identifier,dt=2.05,
                                                     nStates=6,niter=200){
  local_coordination <- rep(NA,length(jobset_str_list))
  for (i in seq_along(jobset_str_list)){
    jobset_str <- jobset_str_list[i]
    job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
    edited_job_id = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
    path_to_est <- file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',
                                                 edited_job_id,'.rds',sep=''))
    if (file.exists(path_to_est)){
      estimate <- readRDS(file=path_to_est)
      local_coordination[i] <- make_local_coordination_agreement_figure_single_cell(jobset_str,estimate,job_id,dt=dt,
                                                     nStates=nStates,niter=niter)
    }
  }
p1 <- ggplot(data=tibble(qq=local_coordination),
       aes(qq)) +
  geom_histogram() +
  theme_bw() +
  labs(x="Quantile in comparison between data and simulations",
       y="Number of cells")
p1
ggsave(here::here("plots/agreement_data_vs_simulations_comparison.eps"),width=6,height=4)
return(local_coordination)
}

make_local_coordination_agreement_figure_single_cell <- function(jobset_str,estimate,job_id,dt=2.05,
                                                     nStates=6,niter=200){

  Data <- process_jobset(jobset_str,K=Inf,max_missing=0.25) %>%
    filter(!is.na(SisterID)) #omit unpaired KTs
  
  sigma_sim <- extract_hidden_states(estimate)

#model was only to fit where sisters had separated by more than 2um
  has_separated_df <- Data %>%
    group_by(SisterPairID,Frame) %>%
    arrange(SisterID,SisterPairID,Frame) %>%
    summarise(kk_dist=-diff(Position_1)) %>%
    group_by(SisterPairID) %>%
    summarise(kk_dist_last = last(kk_dist[!is.na(kk_dist)]),
              separated = kk_dist_last>2.0) #cut off of 2um
  pairIDs <- has_separated_df %>% filter(separated) %>% pull(SisterPairID)
#  pairIDs <- unique(Data$SisterPairID)
  nTracks <- length(pairIDs)
  K <- max(Data$Frame)
cat(paste0("Number of tracks: ",nTracks," Size of MCMC output (single iteration): ",dim(sigma_sim[[1]]),"\n"))  
  t_ana_early <- rstan::extract(estimate,pars="t_ana")$t_ana %>% quantile(.,0.025)
  K_meta <- floor(t_ana_early/dt)
  xi_sampled <- array(0,dim=c(nTracks,K,nStates))
  pp <- array(dim=c(nTracks,K,nStates))
  ff <- array(dim=c(nTracks,K,nStates))
  ss <- array(dim=c(nTracks,K,nStates))
  for (iPair in seq_len(nTracks)){
    for (frame in seq_len(K)){
      for (state in seq_len(nStates)) {
        for (mcmc_iter in seq_len(niter)) {
          xi_sampled[iPair,frame,state] <- xi_sampled[iPair,frame,state] + (sigma_sim[[mcmc_iter]][iPair,frame]==state)
        }
        pp[iPair,frame,state] <- pairIDs[iPair]
        ff[iPair,frame,state] <- frame
        ss[iPair,frame,state] <- state
      }
    }
  }
  xi_sampled <- xi_sampled/niter

  intersister_dist_3d_df <- Data %>%
    group_by(SisterPairID,Frame) %>%
    summarise(intersister_dist_3d = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2))

  intersister_dist_3d_mean_df <- intersister_dist_3d_df %>%
    group_by(SisterPairID) %>%
    summarise(intersister_dist_3d_mean = median(intersister_dist_3d,na.rm=T)) #average in time for each pair

  ############
  #try to use some spatial statistics measures to assess spatial autocorrelation
  library(spdep)
  
  df <- intersister_dist_3d_df %>% 
    inner_join(intersister_dist_3d_mean_df,by="SisterPairID") %>%
    mutate(intersister_dist_3d=intersister_dist_3d-intersister_dist_3d_mean) %>% #subtract mean
    dplyr::select(-intersister_dist_3d_mean) %>%
 #   inner_join(df2 %>% mutate(Frame=frame)) %>% 
    filter(Frame<=K_meta) %>%
    filter(SisterPairID %in% pairIDs) %>% #exclude pairs that have not separated as assessed previously
    # tidyr::gather(stat,value,-SisterPairID,-frame,-Frame) %>%
    inner_join(Data %>% 
                 group_by(SisterPairID,Frame) %>%
                 summarise(x=mean(Position_1),
                           y=mean(Position_2),
                           z=mean(Position_3)))
  
  ###########
  #is a KT pair in the same state as its neighbours
  library(Matrix)
  av_agreement <- rep(0,K_meta)
  all_pairIDs <- has_separated_df %>% filter(separated) %>% pull(SisterPairID)
#  all_pairIDs <- df$SisterPairID %>% unique()
  for (iFrame in seq_len(K_meta)){
    simplified_df <- df %>% filter(Frame==iFrame) %>%
      tidyr::drop_na()
    pairIDs <- unique(simplified_df$SisterPairID)
#    pairIDs <- pairIDs[pairIDs %in% all_pairIDs] #exclude pair if it has not separated as assessed previously
    positions <- cbind(simplified_df[['y']], simplified_df[['z']])
    S_nb <- knn2nb(knearneigh(positions, k=4), row.names=pairIDs)
    # plot(positions)
    # plot(S_nb, positions, add=TRUE, pch=".")
    nb_B <- nb2listw(S_nb, style="B", zero.policy=TRUE)
    sparse_B <- as(nb_B, "CsparseMatrix")
    B <- as.matrix(sparse_B)
    
    for (i in seq_along(pairIDs)){
      reference_kt_state <- xi_sampled[which(pairIDs[i]==all_pairIDs),iFrame,]
      neighbour_indices <- which(B[i,]>0)
      neighbour_state_agreement <- purrr::map_dbl(neighbour_indices,
                                                  function(x) sum(xi_sampled[which(pairIDs[x]==all_pairIDs),iFrame,]*reference_kt_state))
      av_agreement[iFrame] <- av_agreement[iFrame]+mean(neighbour_state_agreement)
    }
    av_agreement[iFrame] <- av_agreement[iFrame]/length(pairIDs)
  }
  
  ######################
  
  
  # I suspect transitions in hidden state can result in the oscillations 
  # ***without any interaction between KT pairs??***
  #  Test this via a simpler simulation model just involving a four state
  #  markov chain
  #  JUH 2021-08-03
  set.seed(123)
  npairs = length(all_pairIDs)
  tLen = 600
  dt = 2
  nstates=4
  
  p_coh = rstan::extract(estimate,pars="p_icoh")$p_icoh %>% mean()
  p_icoh = rstan::extract(estimate,pars="p_coh")$p_coh %>% mean()
  q_coh=1-p_coh
  q_icoh=1-p_icoh
  
  # transition matrix
  M = matrix(c(p_icoh*p_icoh, p_icoh*q_icoh, q_icoh*p_icoh, q_icoh*q_icoh,
               p_coh*q_coh, p_coh*p_coh, q_coh*q_coh, p_coh*q_coh,
               p_coh*q_coh, q_coh*q_coh, p_coh*p_coh, q_coh*p_coh,
               q_icoh*q_icoh, p_icoh*q_icoh, q_icoh*p_icoh, p_icoh*p_icoh),nrow=4,byrow=TRUE)    
  e <- eigen(t(M))
  stationary_distn <- e$vectors[,1]/sum(e$vectors[,1])
  
  nSim <- 100
  ensemble <- rep(NA,nSim)
  av_acf <- array(NA,dim=c(nSim,min(51,K_meta)))
  for (iSim in seq_len(nSim)){
    states = array(0,dim=c(nstates,npairs,tLen/dt+1)) #which state, 46 pairs: ++, +-, -+, --
    # initial states
    which_state=sample.int(nstates,size=46,replace=TRUE,prob=stationary_distn) #randomly sample initial state, equally likely
    for (j in seq_len(npairs)){
      states[which_state[j],j,1] = 1
    }
    
    for (i in seq_len(tLen/dt)){
      for (j in seq_len(npairs)){
        states[,j,i+1] = stats::rmultinom(1, 1, as.numeric(states[,j,i]%*%M)) #mnrnd(1,states(:,j,i)'*M)'
      }
    }
    
    a = apply(states,c(1,3),mean)
    simulation_results_df <- tibble(prob=as.numeric(a),
                                    state=as.numeric(do.call(rbind,
                                                             purrr::map(1:nstates,function(x) rep(x,dim(states)[3])))),
                                    frame=as.numeric(do.call(cbind,
                                                             purrr::map(seq_len(dim(states)[3]),function(x) rep(x,nstates))))
    )
    # simulation_results_df %>%
    # ggplot(aes(frame,prob,color=factor(state))) +
    #   geom_step() +
    #   theme_bw() +
    #   labs(x="Time (s)",
    #        y="Proportion of pairs in each state")
    # # legend({'++','+-','-+','--'});
    
    
    # TODO: next, take locations in metaphase plate from real data, 
    simulated_av_agreement <- rep(0,K_meta)
    all_pairIDs <- df$SisterPairID %>% unique()
    iFrame=1 #get positions and neighbours from first fframe and assume remains similar
    simplified_df <- df %>% filter(Frame==iFrame) %>%
      tidyr::drop_na()
    pairIDs <- unique(simplified_df$SisterPairID)
    positions <- cbind(simplified_df[['y']], simplified_df[['z']])
    S_nb <- knn2nb(knearneigh(positions, k=4), row.names=pairIDs)
    # plot(positions)
    # plot(S_nb, positions, add=TRUE, pch=".")
    nb_B <- nb2listw(S_nb, style="B", zero.policy=TRUE)
    sparse_B <- as(nb_B, "CsparseMatrix")
    B <- as.matrix(sparse_B)
    for (iFrame in seq_len(K_meta)){
      for (i in seq_along(pairIDs)){
        reference_kt_state <- states[,i,iFrame]
        neighbour_indices <- which(B[i,]>0)
        neighbour_state_agreement <- purrr::map_dbl(neighbour_indices,
                                                    function(x) sum(states[,x,iFrame]*reference_kt_state))
        simulated_av_agreement[iFrame] <- simulated_av_agreement[iFrame]+mean(neighbour_state_agreement)
      }
      simulated_av_agreement[iFrame] <- simulated_av_agreement[iFrame]/length(pairIDs)
    }
    # plot(simulated_av_agreement)
    av_acf[iSim,] <- as.numeric(acf(simulated_av_agreement,
                                    lag.max = min(50,K_meta),plot=FALSE)$acf)
    ensemble[iSim] <- mean(simulated_av_agreement)
  }
  ggplot(tibble(autocorr=apply(av_acf,2,mean),
                lag=seq(from=0,to=dt*min(50,K_meta-1),by=dt)),
         aes(x=lag,y=autocorr)) + 
    geom_line() + 
    theme_bw() + labs(x="Lag (s)",
                      y="Autocorrelation")
  ggsave(here::here(paste0("plots/average_autocorrelation_100sims_4state_markov_model_",job_id,".eps")),
         width=4,height=4)
  
  ggplot(tibble(J=ensemble),aes(J)) + 
    geom_histogram() + 
    geom_vline(aes(xintercept=mean(av_agreement)),
               color="red",
               linetype="dashed") + 
    theme_bw() + 
    labs(x="Average alignment in states\nbetween a KT and its neighbours",
         y="Number of simulations")
  ggsave(here::here(paste0("plots/simulated_alignment_vs_observed_alignment_",job_id,".eps")),
         width=210,height=140,units="mm")
  
  pdf(here::here(paste0("plots/positions_of_kt_pairs_in_plate_",job_id,".pdf")))
  plot(positions,xlab='y position (um)',ylab='z position (um)')
  plot(S_nb, positions, add=TRUE, pch=".")
  dev.off()
  
  a <- acf(av_agreement,lag.max = 50,plot=F)
  lag_x = a$lag
  acf_y = a$acf 
  pdf(here::here(paste0("plots/autocorrelation_of_alignment_data_",job_id,".pdf")))
  # acf(av_agreement,lag.max = 50)
  plot(x=lag_x*dt,y=acf_y,type='h',ylab="Autocorrelation",
       xlab="Lag (s)")
  dev.off()
  
  pdf(here::here(paste0("plots/av_agreement_data_time_series_",job_id,".pdf")))
  plot(seq(from=0,to=(length(av_agreement)-1)*dt,by=dt),
       av_agreement,xlab="Time (s)",ylab=expression(paste("Alignment ","J"^"t")))
  lines(seq(from=0,to=(length(av_agreement)-1)*dt,by=dt),
        av_agreement)
  dev.off()
  ######################
  
  # nNeighbours <- 4
  # temp <- array(NA,dim=c(nstates,1+nNeighbours,nSim))
  # out <- rep(NA,nSim)
  # for (iSim in seq_len(nSim)){
  #   for (k in seq_len(1+nNeighbours)) {
  #     temp[,k,iSim] <- stats::rmultinom(1, 1, stationary_distn)
  #   } 
  #   neighbour_state_agreement <- purrr::map_dbl(1:nNeighbours,
  #                                               function(x) sum(temp[,x+1,iSim]*temp[,1,iSim]))
  #   out[iSim] <- mean(neighbour_state_agreement)
  # }
  # hist(out)

#  return(list(mean(av_agreement),av_agreement,ensemble))
sim_ecdf <- stats::ecdf(ensemble)
return(sim_ecdf(mean(av_agreement)))
}
