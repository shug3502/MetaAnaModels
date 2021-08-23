find_switch_pattern <- function(ind,pattern_seq=c(2,1,3),is_joint=FALSE){
  states <- paste(ind,collapse="") #convert into a string
  if (!is_joint){
    stopifnot(length(pattern_seq)==3) #would need to adjust o/w as assume switch preceeded by 2 previous time pts
    loc <- stringr::str_locate_all(states,pattern=paste0(pattern_seq[1],"{2,}",
                                                     pattern_seq[2],"{1,}",
                                                     pattern_seq[3],"{2,}"))
  } else {
    stopifnot(length(pattern_seq)==2)
    loc <- stringr::str_locate_all(states,pattern=paste0(pattern_seq[1],"{2,}",
                                                     pattern_seq[2],"{2,}"))
  } 

  prob <- rep(0,length(ind))
  for (i in seq_along(loc)){
    for (j in seq_len(nrow(loc[[i]]))){
      switch_sequence_str <- paste(ind[loc[[i]][j,1]:loc[[i]][j,2]],collapse="") #just part of sequence that matches the big pattern. ie coherent run, incoherent switch and coherent run
      switch_loc <- stringr::str_locate(switch_sequence_str,pattern=paste(pattern_seq[1:2],collapse="")) #by definition should only have one switch in this
    prob[loc[[i]][j,1] + switch_loc[1]] <- 1 #shift by 1 so mark actual switch point rather than start of pattern
    }
  }
  return(prob)
}

count_num_switch_patterns <- function(ind,pattern_seq=c(2,1,3),is_joint=FALSE){
  states <- paste(ind,collapse="") #convert into a string
  if (!is_joint){
    stopifnot(length(pattern_seq)==3) #would need to adjust o/w as assume switch preceeded by 2 previous time pts
    loc <- stringr::str_locate_all(states,pattern=paste0(pattern_seq[1],"{2,}",
                                                     pattern_seq[2],"{1,}",
                                                     pattern_seq[3],"{2,}"))
  } else {
    stopifnot(length(pattern_seq)==2)
    loc <- stringr::str_locate_all(states,pattern=paste0(pattern_seq[1],"{2,}",
                                                     pattern_seq[2],"{2,}"))
  } 
  num_switches <- nrow(loc[[1]])
  return(num_switches)
}

prob_lead_switch <- function(ind) find_switch_pattern(ind,pattern_seq=c(2,1,3)) +
  find_switch_pattern(ind,pattern_seq=c(3,1,2))
prob_trail_switch <- function(ind) find_switch_pattern(ind,pattern_seq=c(2,4,3)) +
  find_switch_pattern(ind,pattern_seq=c(3,4,2))
prob_joint_switch <- function(ind) find_switch_pattern(ind,pattern_seq=c(2,3),is_joint=TRUE) +
  find_switch_pattern(ind,pattern_seq=c(3,2),is_joint=TRUE)

count_lead_switches <- function(ind) count_num_switch_patterns(ind,pattern_seq=c(2,1,3)) +
  count_num_switch_patterns(ind,pattern_seq=c(3,1,2))
count_trail_switches <- function(ind) count_num_switch_patterns(ind,pattern_seq=c(2,4,3)) +
  count_num_switch_patterns(ind,pattern_seq=c(3,4,2))
count_joint_switches <- function(ind) count_num_switch_patterns(ind,pattern_seq=c(2,3),is_joint=TRUE) +
  count_num_switch_patterns(ind,pattern_seq=c(3,2),is_joint=TRUE)

proportion_lids_estimator <- function(estimate){
  draws <- rstan::extract(estimate,pars="y_tilde",permuted=TRUE)$y_tilde 
  qq <- tibble(lead=apply(draws,1,count_lead_switches),
               trail=apply(draws,1,count_trail_switches),
               joint=apply(draws,1,count_joint_switches))
  qq %>% mutate(total=lead+trail+joint,lids_ratio=lead/total,tids_ratio=trail/total) %>% 
	 summarise(proportion_lids=mean(lids_ratio,na.rm=T),proportion_tids=mean(tids_ratio,na.rm=T),total=mean(total,na.rm=T))
}

plot_hidden_state_draws <- function(draws,nStates=4,dt=2,
                                    possible_states=c("++","+-","-+","--")){
  sz <- dim(draws)
  states <- seq_len(nStates)
  names(states) <- possible_states
  df <- purrr::map_df(states,function(ii) apply(draws,2,function(x) sum(x==ii))/sz[1]) %>%
    mutate(t = seq(from=0,to=(sz[2]-1)*dt,by=dt))
ggplot(df %>% tidyr::gather(state,prob,-t) %>% mutate(state=factor(state,levels=possible_states)), 
       aes(t,prob,fill=factor(state))) + 
  geom_area(color="black",alpha=0.6) + 
  theme_bw()
}

get_switching_df <- function(estimate,K,dt,x){
  draws <- rstan::extract(estimate,pars="y_tilde",permuted=TRUE)$y_tilde  
lead_switches <- apply(draws,1,prob_lead_switch) %>%
  apply(.,1,mean)
trail_switches <- apply(draws,1,prob_trail_switch) %>%
  apply(.,1,mean)
joint_switches <- apply(draws,1,prob_joint_switch) %>%
  apply(.,1,mean)
prop <- proportion_lids_estimator(estimate)
df <- tibble(lead=lead_switches,trail=trail_switches,joint=joint_switches,
             t=seq(from=0,to=(K-1)*dt,by=dt),sis1=x[1,],sis2=x[2,],
             proportion_lids = prop$proportion_lids,
             proportion_tids = prop$proportion_tids,
	     total = prop$total)
return(df)
}

plot_switching_output <- function(estimate,K,dt,x){
  df <- get_switching_df(estimate,K,dt,x)
  draws <- rstan::extract(estimate,pars="y_tilde",permuted=TRUE)$y_tilde  
#visualise
  gg1 <- ggplot(df %>% mutate(total=lead+trail+joint) %>%
                  tidyr::gather(switch,prob,-t,-sis1,-sis2, -proportion_lids),aes(t,prob,color=factor(switch))) +
  geom_line() + theme_bw() +
  theme(legend.position="top") + 
  labs(x="Time (s)",y="Probability\n of switch",title="Probability of LIDS and TIDS switch events",color="Switch type")
gg2 <- ggplot(df %>% tidyr::gather(SisterID,Position,-t,-lead,-trail,-joint,-proportion_lids),
              aes(t,Position,color=factor(SisterID))) + geom_line() +
  theme_bw() +
  theme(legend.position="none") + 
  labs(x="Time (s)",y="Position")

gg3 <- plot_hidden_state_draws(draws,dt=dt) + theme_bw() + 
  labs(x="Time (s)",y="Probability",title="Sampled states",fill="")
#gg4 <- plot_hidden_states(estimate,dt=dt)
print(gg1 + gg2 + gg3 + plot_layout(ncol=1))
}  


get_switching_df_for_track <- function(id,Data,pairIDs,out,K,dt){
  x <- rbind(Data %>% filter(SisterPairID==pairIDs[id]) %>% filter(SisterID==1) %>% pull(Position_1),
             Data %>% filter(SisterPairID==pairIDs[id]) %>% filter(SisterID==2)%>% pull(Position_1))
  #plot_switching_output(out[[id]],K,dt,x)
  df <- get_switching_df(out[[id]],K,dt,x)
  return(df)
}

