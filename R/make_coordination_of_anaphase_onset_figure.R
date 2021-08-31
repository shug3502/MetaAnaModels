make_coordination_of_anaphase_onset_figure <- function(jobset_str,estimate,dt=2.05){

data_single_pair <- process_jobset(jobset_str,K=Inf,max_missing=0.25) %>%
    filter(!is.na(SisterID)) #omit unpaired KTs
  
anaphase_onset_df <- spread_draws(estimate,t_ana[pair]) %>% 
  mutate(SisterPairID = pairIDs[as.integer(pair)]) %>%
  group_by(SisterPairID) %>%
  summarise(t_ana=median(t_ana,na.rm=T))

positions_at_anaphase_onset_df <- anaphase_onset_df %>%
  mutate(Frame=floor(t_ana/dt)) %>%
  left_join(data_single_pair) %>%
  group_by(SisterPairID) %>%
  summarise(t_ana=first(t_ana),
            Frame=first(Frame),
            y=mean(Position_2,na.rm=T),
            z=mean(Position_3,na.rm=T)) %>%
  arrange(t_ana)
  
dist_between_successive_pairs <- sqrt(diff(positions_at_anaphase_onset_df$y)^2 + 
  diff(positions_at_anaphase_onset_df$z)^2)
hist(dist_between_successive_pairs)
av_dist_between_successive_pairs <- median(dist_between_successive_pairs,na.rm=T)


#### randomly permute order of pairs
set.seed(42)
nSim <- 1000
av_dist_between_successive_pairs_sim <- rep(NA,nSim)
for (iSim in seq_len(nSim)){
dummy_df <- positions_at_anaphase_onset_df %>% 
  mutate(dummy_var=rnorm(length(pairIDs))) %>%
  arrange(dummy_var)
dist_between_successive_pairs_sim <- sqrt(diff(dummy_df$y)^2 + 
                                        diff(dummy_df$z)^2)
# hist(dist_between_successive_pairs_sim)
av_dist_between_successive_pairs_sim[iSim] <- median(dist_between_successive_pairs_sim,na.rm=T)
}
ggplot(data=tibble(dist=av_dist_between_successive_pairs_sim),
       aes(dist)) + 
  geom_histogram() +
  geom_vline(xintercept=av_dist_between_successive_pairs,color="red",linetype="dashed") + 
  theme_bw() + 
  labs(x="Average distance between kinetochore\npairs at anaphase onset (um)",
       y="Number of simulations")
ggsave(here::here("plots/av_distance_between_pairs_at_anaphase_onset.eps"),
       width=6,height=4)
}
