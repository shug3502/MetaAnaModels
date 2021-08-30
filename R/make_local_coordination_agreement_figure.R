make_local_coordination_agreement_figure <- function(jobset_str,estimate,dt=2.05,
                                                     nStates=6,niter=200){

  data_single_pair <- process_jobset(jobset_str,K=Inf,max_missing=0.25) %>%
    filter(!is.na(SisterID)) #omit unpaired KTs
  
  sigma_sim <- extract_hidden_states(estimate)

  pairIDs <- unique(data_single_pair$SisterPairID)
  nTracks <- length(pairIDs)
  K <- max(data_single_pair$Frame)
  
  t_ana_early <- rstan::extract(estimate,pars="t_ana")$t_ana %>% quantile(.,0.025)
  K_meta <- floor(t_ana_early/dt)
  # xi_sampled <- array(0,dim=c(nTracks,K,nStates))
  # pp <- array(dim=c(nTracks,K,nStates))
  # ff <- array(dim=c(nTracks,K,nStates))
  # ss <- array(dim=c(nTracks,K,nStates))
  # for (iPair in seq_len(nTracks)){
  #   for (frame in seq_len(K)){
  #     for (state in seq_len(nStates)) {
  #       for (mcmc_iter in seq_len(niter)) {
  #         xi_sampled[iPair,frame,state] <- xi_sampled[iPair,frame,state] + (sigma_sim[[mcmc_iter]][iPair,frame]==state)
  #       }
  #       pp[iPair,frame,state] <- pairIDs[iPair]
  #       ff[iPair,frame,state] <- frame
  #       ss[iPair,frame,state] <- state
  #     }
  #   }
  # }
  # xi_sampled <- xi_sampled/niter
  # df <- tibble(xi=as.numeric(xi_sampled),
  #              SisterPairID=as.integer(pp),
  #              frame=as.integer(ff),state=as.integer(ss))
  # df %>% ggplot(aes(frame,xi,color=factor(state))) + geom_step() +
  #   facet_wrap(.~SisterPairID,labeller = label_both) + theme_bw()
  # 
  ################

  intersister_dist_3d_df <- data_single_pair %>%
    group_by(SisterPairID,Frame) %>%
    summarise(intersister_dist_3d = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2))

  intersister_dist_3d_mean_df <- intersister_dist_3d_df %>%
    group_by(SisterPairID) %>%
    summarise(intersister_dist_3d_mean = median(intersister_dist_3d,na.rm=T)) #average in time for each pair

  # ggplot(intersister_dist_3d_df %>% filter(Frame<240),
  #        aes((Frame-1)*dt,intersister_dist_3d,group=SisterPairID)) + 
  #   geom_line() + theme_bw() + 
  #   facet_wrap(.~ SisterPairID)
  
  df2 <- tibble(xi=as.numeric(xi_sampled[,,2]-xi_sampled[,,3]),
                SisterPairID=as.integer(pp[,,2]),
                frame=as.integer(ff[,,2]))
  # 
  # intersister_dist_3d_df %>% 
  #   inner_join(intersister_dist_3d_mean_df,by="SisterPairID") %>%
  #   mutate(intersister_dist_3d=intersister_dist_3d-intersister_dist_3d_mean) %>% #subtract mean
  #   dplyr::select(-intersister_dist_3d_mean) %>%
  #   inner_join(df2 %>% mutate(Frame=frame)) %>% filter(Frame<240) %>%
  #   tidyr::gather(stat,value,-SisterPairID,-frame,-Frame) %>%
  #   ggplot(aes(Frame,value,color=stat)) + 
  #   geom_line() +
  #   facet_wrap(.~SisterPairID,labeller = label_both) + theme_bw()
  
  # 
  # library(gganimate)
  # plt <- data_single_pair %>% group_by(SisterPairID) %>%
  #   summarise(av_y=median(Position_2,na.rm=T),
  #             av_z=median(Position_3,na.rm=T)) %>%
  #   mutate(y=cut(av_y,4),
  #          z=cut(av_z,4)) %>%
  #   inner_join(intersister_dist_3d_df %>% 
  #                inner_join(intersister_dist_3d_mean_df,by="SisterPairID") %>%
  #                mutate(intersister_dist_3d=intersister_dist_3d-intersister_dist_3d_mean) %>% #subtract mean
  #                dplyr::select(-intersister_dist_3d_mean) %>%
  #                inner_join(df2 %>% mutate(Frame=frame)) %>% filter(Frame<240)) %>%
  #   #tidyr::gather(stat,value,-SisterPairID,-frame,-Frame) %>%
  #   ggplot(aes(av_y,av_z,color=xi)) + 
  #   geom_point() +
  #   theme_bw()
  # 
  # anim <- plt + transition_states(Frame,
  #                                 transition_length = 0.1,
  #                                 state_length = 1) +
  #   ggtitle('Frame {frame} of {nframes}')
  # 
  # animate(anim,nframes=K)
  # anim_save(here::here("plots/xi_animation.gif"))
  
  ############
  #try to use some spatial statistics measures to assess spatial autocorrelation
  library(spdep)
  
  df <- intersister_dist_3d_df %>% 
    inner_join(intersister_dist_3d_mean_df,by="SisterPairID") %>%
    mutate(intersister_dist_3d=intersister_dist_3d-intersister_dist_3d_mean) %>% #subtract mean
    dplyr::select(-intersister_dist_3d_mean) %>%
    inner_join(df2 %>% mutate(Frame=frame)) %>% 
    filter(Frame<240) %>%
    # tidyr::gather(stat,value,-SisterPairID,-frame,-Frame) %>%
    inner_join(data_single_pair %>% 
                 group_by(SisterPairID,Frame) %>%
                 summarise(x=mean(Position_1),
                           y=mean(Position_2),
                           z=mean(Position_3)))
  
  # get_moran_I <- function(frame_num,within_dist=1.5){
  #   simplified_df <- df %>% filter(Frame==frame_num) %>%
  #   tidyr::drop_na()
  # # kk_dist_vec <- simplified_df$intersister_dist_3d
  # positions <- cbind(simplified_df[['y']], simplified_df[['z']])
  # pairwise_dist <- dist(positions) %>% as.matrix()
  # weights <- 1/pairwise_dist
  # diag(weights) <- 0
  # 
  # # ape::Moran.I(simplified_df$xi, weights)
  # if (length(simplified_df$xi)==0) return(NA)
  # moran_i <- ape::Moran.I(simplified_df$xi, (pairwise_dist > 0 & (pairwise_dist <= within_dist)))
  # # ape::Moran.I(simplified_df$xi, (pairwise_dist > 0 & pairwise_dist <= 1.5))
  # return(moran_i$observed)
  # }
  # out <- purrr::map_dbl(1:239,get_moran_I)
  # plot(out)
  # 
  # get_correlogram <- function(frame_num){
  #   simplified_df <- df %>% filter(Frame==frame_num) %>%
  #     tidyr::drop_na()
  #   positions <- cbind(simplified_df[['y']], simplified_df[['z']])
  #  
  #   # S_nb <- knn2nb(knearneigh(positions, k=5), row.names=unique(simplified_df$SisterPairID))
  #   # dsts <- unlist(nbdists(S_nb, positions))
  #   # max_1nn <- max(dsts)
  #   max_1nn <- 3
  #   S_nb <- dnearneigh(positions, d1=0, d2=1*max_1nn, row.names=unique(simplified_df$SisterPairID))
  #   correl <- sp.correlogram(S_nb, simplified_df$xi,order = 2)
  #   return(correl$res)
  # }
  # correlogram.out <- purrr::map(1:239,get_correlogram) 
  # purrr::map_dbl(1:2,function(y) mean(purrr::map_dbl(1:239,function(x) correlogram.out[[x]][y])))
  # plot(out)
  # 
  # S_nb <- knn2nb(knearneigh(positions, k=5), row.names=unique(simplified_df$SisterPairID))
  # plot(positions)
  # plot(S_nb, positions, add=TRUE, pch=".")
  # sp.correlogram(S_nb, simplified_df$xi,order = 3)
  # moran.plot(simplified_df$xi, nb2listw(S_nb),label=F)
  # 
  # dsts <- unlist(nbdists(S_nb, positions))
  # max_1nn <- max(dsts)
  # S11_nb <- dnearneigh(positions, d1=0, d2=0.75*max_1nn, row.names=unique(simplified_df$SisterPairID))
  # S12_nb <- dnearneigh(positions, d1=0, d2=1*max_1nn, row.names=unique(simplified_df$SisterPairID))
  # plot(positions)
  # plot(S11_nb, positions, add=TRUE, pch=".")
  # 
  # plot(positions)
  # plot(S12_nb, positions, add=TRUE, pch=".")
  # 
  # sp.correlogram(S11_nb, simplified_df$xi,order = 3)
  # moran.plot(simplified_df$xi, nb2listw(S11_nb),label=F)
  
  
  #############
  
  pairwise_correlation_KT_pairs <- function(df,pair_i,pair_j){
    simplified_df_i <- df %>% filter(SisterPairID==pair_i) %>% drop_na()
    simplified_df_j <- df %>% filter(SisterPairID==pair_j) %>% drop_na()
    mask <- intersect(simplified_df_i$Frame,simplified_df_j$Frame)
    xx <- simplified_df_i %>% filter(Frame %in% mask) %>% pull(xi)
    yy <- simplified_df_j %>% filter(Frame %in% mask) %>% pull(xi)
    # positions <- cbind(simplified_df[['y']], simplified_df[['z']])
    cor(xx,yy)
  }
  
  pairwise_dist_between_KT_pairs <- function(df,pair_i,pair_j){
    av_pos_i <- df %>% filter(SisterPairID==pair_i) %>% 
      summarise(y=median(y,na.rm=T),
                z=median(z,na.rm=T))
    av_pos_j <- df %>% filter(SisterPairID==pair_j) %>% 
      summarise(y=median(y,na.rm=T),
                z=median(z,na.rm=T))  
    sqrt((av_pos_i$y-av_pos_j$y)^2 + (av_pos_i$z-av_pos_j$z)^2)
  }
  
  all_pairs <- expand.grid(i=pairIDs,j=pairIDs)
  all_cors <- purrr::map2_dbl(all_pairs$i,all_pairs$j,function(x1,y1) pairwise_correlation_KT_pairs(df,x1,y1))
  all_dists <- purrr::map2_dbl(all_pairs$i,all_pairs$j,function(x2,y2) pairwise_dist_between_KT_pairs(df,x2,y2))
  ggplot(tibble(dist=all_dists,cors=all_cors),aes(dist,cors)) +
    geom_point() + geom_smooth() + theme_bw()
  
  # av_pos_df <- df %>% group_by(SisterPairID) %>%
  #   summarise(av_y=median(y,na.rm=T),av_z=median(z,na.rm=T))
  # df_with_dist_from_ref_pair_df <- av_pos_df %>% 
  #   mutate(d=sqrt((av_y-last(av_pos_df$av_y))^2+(av_z-last(av_pos_df$av_z))^2)) %>%
  #   mutate(discrete_d=cut(d,5)) %>%
  #   inner_join(df)
  # ggplot(df_with_dist_from_ref_pair_df,
  #        aes((Frame-1)*dt,xi,color=d,group=SisterPairID)) + 
  #   geom_line() + 
  #   geom_line(data=df_with_dist_from_ref_pair_df %>% filter(SisterPairID==last(pairIDs)),
  #             aes((Frame-1)*dt,xi),color="green") +
  #   theme_bw() + facet_wrap(.~discrete_d)
  
  
  ###########
  #is a KT pair in the same state as its neighbours
  library(Matrix)
  av_agreement <- rep(0,K_meta)
  all_pairIDs <- df$SisterPairID %>% unique()
  for (iFrame in seq_len(K_meta)){
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
  av_acf <- array(NA,dim=c(nSim,51))
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
                                    lag.max = 50,plot=FALSE)$acf)
    ensemble[iSim] <- mean(simulated_av_agreement)
  }
  ggplot(tibble(autocorr=apply(av_acf,2,mean),
                lag=seq(from=0,to=100,by=dt)),
         aes(x=lag,y=autocorr)) + 
    geom_line() + 
    theme_bw() + labs(x="Lag (s)",
                      y="Autocorrelation")
  ggsave(here::here("plots/average_autocorrelation_100sims_4state_markov_model.eps"),device=cairo_ps,
         width=4,height=4)
  
  ggplot(tibble(J=ensemble),aes(J)) + 
    geom_histogram() + 
    geom_vline(aes(xintercept=mean(av_agreement)),
               color="red",
               linetype="dashed") + 
    theme_bw() + 
    labs(x="Average agreement in states\nbetween a KT and its neighbours",
         y="Number of simulations")
  ggsave(here::here("plots/simulated_agreement_vs_observed_agreement.eps"),
         device=cairo_ps,width=210,height=140,units="mm")
  
  pdf(here::here("plots/positions_of_kt_pairs_in_plate.pdf"))
  plot(positions,xlab='y position (um)',ylab='z position (um)')
  plot(S_nb, positions, add=TRUE, pch=".")
  dev.off()
  
  a <- acf(av_agreement,lag.max = 50,plot=F)
  lag_x = a$lag
  acf_y = a$acf 
  pdf(here::here("plots/autocorrelation_of_alignment_data.pdf"))
  # acf(av_agreement,lag.max = 50)
  plot(x=lag_x*dt,y=acf_y,type='h',ylab="Autocorrelation",
       xlab="Lag (s)")
  dev.off()
  
  pdf(here::here("plots/av_agreement_data_time_series.pdf"))
  plot(seq(from=0,to=(length(av_agreement)-1)*dt,by=dt),
       av_agreement,xlab="Time (s)",ylab=expression(paste("Agreement ","J"^"t")))
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
  
}