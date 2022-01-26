make_coordination_of_anaphase_onset_figure <- function(identifier,jobset_str_list,dt=2.05){
empirical_quantile_of_obs_vs_sim <- rep(NA,length(jobset_str_list))
#av_diff_from_hist <- rep(0,30)
nPairs_at_ana_vec <- rep(NA,length(jobset_str_list))
nPts_at_ana <- rep(NA,length(jobset_str_list))
pdist_t_ana_df <- data.frame()
  for (i in seq_along(jobset_str_list)){
print(i)
print(jobset_str_list[i])
    job_id = stringr::str_split(jobset_str_list[i],"kittracking")[[1]][2]
    edited_job_id = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
    path_to_est <- file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',
                                                 edited_job_id,'.rds',sep=''))
    if (file.exists(path_to_est)){
      estimate <- readRDS(file=path_to_est)
      lst <- make_coordination_of_anaphase_onset_figure_single_cell(jobset_str_list[i],estimate,edited_job_id,dt=dt)
      empirical_quantile_of_obs_vs_sim[i] <- lst[[1]]
      nPts_at_ana[i] <- lst[[2]]
#      av_diff_from_hist <- av_diff_from_hist + lst[[2]]
#      pdist_t_ana_df <- rbind(pdist_t_ana_df,lst[[2]])
#      nPairs_at_ana_vec[i] <- lst[[3]]
    } #otherwise skip cell
  }
p1 <- ggplot(data=tibble(dist=empirical_quantile_of_obs_vs_sim),
       aes(dist)) + 
  geom_histogram() +
  theme_bw() +
  labs(x="Percentile in comparison\nbetween data and simulations",
       y="Number of cells")
#p2 <- ggplot(data=tibble(dist=empirical_quantile_of_obs_vs_sim,nPairs=nPairs_at_ana_vec),
#	     aes(nPairs,dist)) + 
#	geom_point() + 
#        theme_bw()
#p1 | p2
p1
ggsave(here::here("plots/av_distance_between_pairs_at_anaphase_onset_quantile_histogram.eps"),width=6,height=4)
p2 <- ggplot(data=tibble(quantile=empirical_quantile_of_obs_vs_sim,
			 nPts=nPts_at_ana),
	     aes(nPts,quantile)) +
	labs(x="Number of tracked pairs at anaphase onset",
	     y="Percentile in comparison between data and simulations") +  
	geom_vline(xintercept=30,linetype="dashed",color="grey") +
	geom_point() + theme_bw()
p2
ggsave(here::here("plots/number_of_points_in_square_vs_quantile.eps"),width=6,height=4)
p3 <- ggplot(data=tibble(quantile=empirical_quantile_of_obs_vs_sim,
                         nPts=nPts_at_ana) %>% filter(nPts>30),
       aes(quantile)) +
  geom_histogram() +
  theme_bw() +
  labs(x="Percentile in comparison\nbetween data and simulations",
       y="Number of cells")
p2 / p3
ggsave(here::here("plots/number_of_points_in_square_vs_quantile_plus_restricted_to_more_than_30_pairs.eps"),width=6,height=4)

#plt <- ggplot(pdist_t_ana_df,aes(dist,t_ana_diff)) + 
#geom_point(size=0.5) + 
#geom_smooth(color="black") +
#theme_bw()
#plt
#ggsave(here::here("plots/spatial_all_the_cells.eps"),width=6,height=4)

#av_diff_from_hist <- av_diff_from_hist/length(jobset_str_list)
#ggplot(data=data.frame(x=seq(from=0.25,to=14.75,by=0.5),y=av_diff_from_hist),
#       aes(x,y)) + 
#geom_bar(stat="identity") + 
#theme_bw()
#ggsave(here::here("plots/av_distance_between_pairs_distance_between_histograms.eps"),width=6,height=4)
return(empirical_quantile_of_obs_vs_sim)
}

#my_kernel <- function(delt,dely,delz,lengthscale,timescale){
#(1/(2*pi)*exp(-.5 * (dely^2+delz^2)/lengthscale^2)/lengthscale^2)*(timescale*exp(-timescale*abs(delt)))
#}

#compute_av_kernel_on_all_interactions <- function(positions_at_anaphase_onset_df,method="random",
#						  lengthscale=2,timescale=5){
#  if (method=="random"){
#    #permute times of anaphase onset
#    positions_at_anaphase_onset_df <- positions_at_anaphase_onset_df %>% mutate(t_ana = sample(t_ana,length(t_ana),replace=FALSE))
#  }
#  nSisters <- nrow(positions_at_anaphase_onset_df)
#  k <- matrix(0,nrow=nSisters,ncol=nSisters)
#  for (i in 1:nSisters){
#    for (j in setdiff(1:nSisters,i)){
#      k[i,j] <- my_kernel(positions_at_anaphase_onset_df[i,"t_ana"][[1]]-positions_at_anaphase_onset_df[j,"t_ana"][[1]],
#			positions_at_anaphase_onset_df[i,"y"][[1]]-positions_at_anaphase_onset_df[j,"y"][[1]],
#			positions_at_anaphase_onset_df[i,"z"][[1]]-positions_at_anaphase_onset_df[j,"z"][[1]],
#			lengthscale,timescale)
#    }
#  }
#  av_k <- sum(k)/(nSisters*(nSisters-1))
#}

make_coordination_of_anaphase_onset_figure_single_cell <- function(jobset_str,estimate,job_id,dt=2.05){

Data <- process_jobset(jobset_str,K=Inf,max_missing=0.25) %>%
    filter(!is.na(SisterID)) #omit unpaired KTs

#model was only to fit where sisters had separated by more than 2um
  has_separated_df <- Data %>%
    group_by(SisterPairID,Frame) %>%
    arrange(SisterID,SisterPairID,Frame) %>%
    summarise(kk_dist=-diff(Position_1)) %>%
    group_by(SisterPairID) %>%
    summarise(kk_dist_last = last(kk_dist[!is.na(kk_dist)]),
              separated = kk_dist_last>2.0) #cut off of 2um
pairIDs <- has_separated_df %>% filter(separated) %>% pull(SisterPairID)
  
anaphase_onset_df <- spread_draws(estimate,t_ana[pair]) %>% 
  mutate(SisterPairID = pairIDs[as.integer(pair)]) %>%
  group_by(SisterPairID) %>%
  summarise(t_ana=median(t_ana,na.rm=T))

positions_at_anaphase_onset_df <- anaphase_onset_df %>%
  mutate(Frame=floor(t_ana/dt)) %>%
  left_join(Data) %>%
  group_by(SisterPairID) %>%
  summarise(t_ana=first(t_ana),
            Frame=first(Frame),
            y=mean(Position_2,na.rm=T),
            z=mean(Position_3,na.rm=T)) %>%
  arrange(t_ana)
  plt1 <- ggplot(positions_at_anaphase_onset_df,
		 aes(y,z,color=t_ana)) + 
	  geom_point() + geom_path() + theme_bw()

empirical_quantile <- ana_spatial_analysis(positions_at_anaphase_onset_df,job_id,nRepeats=1000)

#pdist_t_ana_df <- get_spatial_corr(positions_at_anaphase_onset_df,job_id)
#dist_between_successive_pairs <- sqrt(diff(positions_at_anaphase_onset_df$y)^2 + 
#  diff(positions_at_anaphase_onset_df$z)^2)
#hist(dist_between_successive_pairs)
#av_dist_between_successive_pairs <- median(dist_between_successive_pairs,na.rm=T)
#dist_between_pairs_two_apart <- sqrt((positions_at_anaphase_onset_df$y - lag(positions_at_anaphase_onset_df$y,2))^2 +
#  (positions_at_anaphase_onset_df$z-lag(positions_at_anaphase_onset_df$z,2))^2)
#av_dist_between_pairs_two_apart <- median(dist_between_pairs_two_apart,na.rm=T)

#### randomly permute order of pairs
#set.seed(42)
#nPairs_at_ana <- length(unique(positions_at_anaphase_onset_df$SisterPairID))
#stopifnot(nPairs_at_ana>10)
#nSim <- 1000
#av_dist_between_successive_pairs_sim <- rep(NA,nSim)
#all_sim_dists_store <- c()
#for (iSim in seq_len(nSim)){
#dummy_df <- positions_at_anaphase_onset_df %>% 
#  mutate(dummy_var=rnorm(nPairs_at_ana)) %>%
#  arrange(dummy_var)
#dist_between_successive_pairs_sim <- sqrt(diff(dummy_df$y)^2 + 
#                                        diff(dummy_df$z)^2)
## hist(dist_between_successive_pairs_sim)
#av_dist_between_successive_pairs_sim[iSim] <- median(dist_between_successive_pairs_sim,na.rm=T)
#all_sim_dists_store <- c(all_sim_dists_store,dist_between_successive_pairs_sim)
#}
#if (sum(is.na(av_dist_between_successive_pairs_sim))>nSim/2) {return(NA)}
#plt2 <- ggplot(data=tibble(dist=av_dist_between_successive_pairs_sim),
#       aes(dist)) + 
#  geom_histogram() +
#geom_histogram(data=tibble(dist=dist_between_successive_pairs),aes(y=..count../sum(..count..)),fill="cyan") +
#geom_histogram(data=tibble(dist=all_sim_dists_store),aes(y=..count../sum(..count..)),fill="magenta") +
#  geom_vline(xintercept=av_dist_between_successive_pairs,color="red",linetype="dashed") + 
##  geom_vline(xintercept=av_dist_between_pairs_two_apart,color="magenta",linetype="dashed") +
#  theme_bw() + 
#  labs(x="Average distance between kinetochore\npairs at anaphase onset (um)",
#       y="Number of simulations")
#plt3 <- ggplot(data=tibble(dist=av_dist_between_successive_pairs_sim),
#       aes(dist)) +
#  geom_histogram(aes(y=..count../sum(..count..))) +
#  geom_histogram(data=tibble(dist=dist_between_successive_pairs),aes(y=..count../sum(..count..)),fill="cyan") +
#  geom_histogram(data=tibble(dist=all_sim_dists_store),aes(y=..count../sum(..count..)),fill="magenta") +
#  geom_vline(xintercept=av_dist_between_successive_pairs,color="red",linetype="dashed") +
##  geom_vline(xintercept=av_dist_between_pairs_two_apart,color="magenta",linetype="dashed") +
#  theme_bw() +
#  labs(x="Average distance between kinetochore\npairs at anaphase onset (um)",
#       y="Number of simulations")
#plt1 | plt3
#ggsave(here::here(paste0("plots/av_distance_between_pairs_at_anaphase_onset_",job_id,".eps")),
#       width=6,height=4)
#sim_ecdf <- stats::ecdf(av_dist_between_successive_pairs_sim)

#df = data.frame(y = c(dist_between_successive_pairs, all_sim_dists_store),
#                x = c(rep(0,length(dist_between_successive_pairs)),rep(1,length(all_sim_dists_store))))
##full hist
#fullhist = hist(df$y, breaks = seq(from=0,to=15,by=0.5)) #specify more breaks than probably necessary
##create histograms for 0 & 1 using breaks from full histogram
#zerohist = with(subset(df, x == 0), hist(y, breaks = fullhist$breaks))
#oneshist = with(subset(df, x == 1), hist(y, breaks = fullhist$breaks))
##combine the hists
#combhist = fullhist
#combhist$counts = zerohist$counts/sum(zerohist$counts) - oneshist$counts/sum(oneshist$counts)

#av_kern_data <- compute_av_kernel_on_all_interactions(positions_at_anaphase_onset_df,method="data")
#av_kern_sim <- rep(NA,nSim)
#for (ii in 1:nSim){
#  av_kern_sim[ii] <- compute_av_kernel_on_all_interactions(positions_at_anaphase_onset_df,method="random")
#}
#plt4 <- ggplot(data=tibble(dist=av_kern_sim),
#       aes(dist)) +
#  geom_histogram() +
#  geom_vline(xintercept=av_kern_data,color="red",linetype="dashed") +
#  theme_bw() +
#  labs(x="Kernel based on distance and time between kinetochore\npairs at anaphase onset",
#       y="Number of simulations")
#ggsave(here::here(paste0("plots/kernel_histogram_",job_id,".eps")),
#       width=6,height=4)

#return(sim_ecdf(av_dist_between_successive_pairs))
#return(sim_ecdf(av_dist_between_pairs_two_apart))
#return(list(sim_ecdf(av_dist_between_successive_pairs),combhist$counts,nPairs_at_ana))
#return(list(sim_ecdf(av_dist_between_successive_pairs),pdist_t_ana_df,nPairs_at_ana))
#return(list(empirical_quantile,pdist_t_ana_df,nPairs_at_ana))
return(empirical_quantile)
}
