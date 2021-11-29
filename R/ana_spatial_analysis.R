ana_spatial_analysis <- function(positions_at_anaphase_onset_df,identifier,nRepeats=1000){

nPairs <- nrow(positions_at_anaphase_onset_df)
#from the data
S <- matrix(rep(NA,nPairs*(nPairs-1)),ncol=2)
k <- 0; 
for (i in 1:nPairs){
#  for (j in setdiff(1:nPairs,i)){
  for (j in setdiff(seq_len(i),i)){
    k = k+1;
    S[k,1] <- abs(positions_at_anaphase_onset_df[["t_ana"]][i]-positions_at_anaphase_onset_df[["t_ana"]][j])
    S[k,2] <- sqrt((positions_at_anaphase_onset_df[["y"]][i]-positions_at_anaphase_onset_df[["y"]][j])^2 + (positions_at_anaphase_onset_df[["z"]][i]-positions_at_anaphase_onset_df[["z"]][j])^2)
  }
}

num_bottom_left_resample <- purrr::map_int(1:nRepeats,function(x) resample_distances(positions_at_anaphase_onset_df) %>% get_number_in_bottom_left(.,a=4,b=2))
num_bottom_left_data <- get_number_in_bottom_left(S,a=4,b=2)
sim_ecdf <- stats::ecdf(num_bottom_left_resample)

p1 <- ggplot(data=tibble(tt=S[,1],dd=S[,2]),
	     aes(tt,dd)) + 
	geom_point() + theme_bw() + 
        geom_vline(xintercept=2,linetype="dashed",color="red") +
        geom_hline(yintercept=2,linetype="dashed",color="red") +
labs(x="Difference in anaphase onset time (s)",
     y="Distance between pairs (um)")
S_random_plot <- resample_distances(positions_at_anaphase_onset_df)
q1 <-  ggplot(data=tibble(tt=S_random_plot[,1],dd=S_random_plot[,2]),
             aes(tt,dd)) + 
        geom_point() + theme_bw() +
        geom_vline(xintercept=2,linetype="dashed",color="red") +
        geom_hline(yintercept=2,linetype="dashed",color="red") +
labs(x="Difference in anaphase onset time (s)",
     y="Distance between pairs (um)")
p1 | q1
ggsave(here::here(paste0("plots/dist_vs_time_difference_of_pairs_at_anaphase_",identifier,".eps")))

pdf(here::here(paste0("plots/ecdf_",identifier,".pdf")))
plot(sim_ecdf,ylab="eCDF",xlab="Number of points",main="")
abline(v=num_bottom_left_data,col = "gray60",lty="dashed")
dev.off()

return(list(sim_ecdf(num_bottom_left_data),nPairs,num_bottom_left_data))
}

get_number_in_bottom_left <- function(S,a=2,b=1){
sum((S[,1]<a)&(S[,2]<b))
}

#resample_distances <- function(S){
##resample the distances for comparison
#  S_random <- matrix(rep(NA,2*nrow(S)),ncol=2)
#  S_random[,1] <- S[,1]
#  S_random[,2] <- sample(S[,2],size=nrow(S),replace=FALSE)
#  return(S_random)
#}
resample_distances <- function(positions_at_anaphase_onset_df){
nPairs <- nrow(positions_at_anaphase_onset_df)
positions_at_anaphase_onset_df <- positions_at_anaphase_onset_df %>%
  mutate(t_ana = sample(t_ana,size=nPairs,replace=FALSE))
#from the data
S <- matrix(rep(NA,nPairs*(nPairs-1)),ncol=2)
k <- 0; 
for (i in 1:nPairs){
#  for (j in setdiff(1:nPairs,i)){
  for (j in setdiff(seq_len(i),i)){
    k = k+1;
    S[k,1] <- abs(positions_at_anaphase_onset_df[["t_ana"]][i]-positions_at_anaphase_onset_df[["t_ana"]][j])
    S[k,2] <- sqrt((positions_at_anaphase_onset_df[["y"]][i]-positions_at_anaphase_onset_df[["y"]][j])^2 + (positions_at_anaphase_onset_df[["z"]][i]-positions_at_anaphase_onset_df[["z"]][j])^2)
  }
}
return(S)
}
