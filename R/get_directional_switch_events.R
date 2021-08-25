myf <- function(x,prob_fun,pairIDs,K){ 
  a <- apply(x,1,prob_fun)
  colnames(a) <- pairIDs
  df <- as.data.frame(cbind(a,frame=seq_len(K)))
  return(tidyr::gather(df,SisterPairID,prob,-frame))
}

get_switch_probs_df <- function(sigma_sim,pairIDs,K){
  trail_switch_probs <- purrr::map_df(sigma_sim,function(x) myf(x,prob_trail_switch,pairIDs,K),.id=".iter") %>%
    mutate(switch_type="trail")
  lead_switch_probs <- purrr::map_df(sigma_sim,function(x) myf(x,prob_lead_switch,pairIDs,K),.id=".iter") %>%
    mutate(switch_type="lead")
  joint_switch_probs <- purrr::map_df(sigma_sim,function(x) myf(x,prob_joint_switch,pairIDs,K),.id=".iter") %>%
    mutate(switch_type="joint")
  switch_probs_df <- bind_rows(trail_switch_probs,
          lead_switch_probs,
          joint_switch_probs) %>% as_tibble()
  return(switch_probs_df)
}

get_directional_switch_events <- function(switch_probs_df, Data, dt=2.05, min_prob=0.5){
xyt_df = switch_probs_df %>% 
  group_by(frame,SisterPairID,switch_type) %>%
  summarise(prob=mean(prob)) %>% #average over mcmc iterations
  ungroup() %>%
  filter(prob>min_prob) %>% 
  mutate(Frame=frame,SisterPairID=as.integer(SisterPairID)) %>% 
  inner_join(Data,by=c("SisterPairID","Frame")) %>%
  mutate(Time=(Frame-1)*dt) %>%
  group_by(Frame,SisterPairID) %>%
  summarise(Position_1=mean(Position_1),
            Position_2=mean(Position_2),
            Position_3=mean(Position_3),
            Time=first(Time))
  return(xyt_df)
}
