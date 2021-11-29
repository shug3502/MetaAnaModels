how_common_are_reversals <- function(jobset_str,sigma_sim,niter=2000){
total_reversals <- purrr::map(1:niter,function(iter) apply(sigma_sim[[iter]],1,function(x) sum(x==6))) %>% 
purrr::reduce(`+`)
av_frames_of_reversal <- total_reversals/niter

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
if (length(pairIDs)!=length(av_frames_of_reversal)){ pairIDs <- has_separated_df$SisterPairID} #fall back on using all pairs if not going to work
return(tibble(filename=rep(jobset_str,length(av_frames_of_reversal)),SisterPairID=pairIDs,reversals=av_frames_of_reversal))
}

get_av_radius_of_pairs <- function(jobset_str){
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

Data %>% filter(SisterPairID %in% pairIDs) %>%
group_by(SisterPairID) %>%
summarise(av_radius=median(sqrt(Position_2^2+Position_3^2),na.rm=T)) %>%
mutate(filename=jobset_str)
}

get_laziness_for_cell <- function(jobset_str,draws,min_num_sisters=10,nframes_early=300){
library(lazychromosomes) #https://github.com/shug3502/lazychromosomes
  Data <- process_jobset(jobset_str,K=Inf,max_missing=0.25) %>%
    filter(!is.na(SisterID)) %>% #omit unpaired KTs
    mutate(filename=jobset_str) %>%
    add_kittracking_column()
  start_of_anaphase_df <- draws %>%
    mutate(t_ana = if_else(param=="t_ana",theta,NA_real_)) %>%
    group_by(cell,filename) %>%
    summarise(ana_start = quantile(t_ana,probs=0.5,na.rm=TRUE)) %>% #find time of first separating sister for start of anaphase
    ungroup() %>%
    mutate(cell=as.integer(cell),
           kittracking_file_str = filename %>% stringr::str_split('/') %>% purrr::map(last)) %>%
    tidyr::unnest()
  positions_by_frame_df <- Data %>%
    inner_join(start_of_anaphase_df,by="kittracking_file_str") %>%
    mutate(Frames_since_start = Frame - ana_start%/%dt) %>%
    filter(Frames_since_start < nframes_early) %>%
    filter(Frames_since_start >= -nframes_early) %>%
    group_by(kittracking_file_str,Frames_since_start,SisterPairID,SisterID) %>%
    summarise(Position_1 = first(Position_1),
              Position_2 = first(Position_2),
              Position_3 = first(Position_3)) %>%
    group_by(kittracking_file_str) %>%
    filter(length(unique(SisterPairID))>min_num_sisters)
  positions_by_frame_df %>% get_laziness() %>%
    group_by(SisterPairID) %>%
    summarise(max_laziness=max(laziness,na.rm=T)) %>%
    mutate(filename=jobset_str)
}

