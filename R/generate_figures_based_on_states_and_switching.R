#TODO: loop over jobset_strs in a jobset_str_list. Output some data frames with results that can be combined and plotted together?
get_all_interesting_stats <- function(Data, draws, ana_margin=60, 
                                      window_size=20,dt=2.05, 
                                      t_ana_quantile=0.5,
                                      treatment_str_vec="Unperturbed"){
  #function to extract all the appropriate summary stats
  
  start_of_anaphase_df <- draws %>%
    mutate(t_ana = if_else(param=="t_ana",theta,NA_real_)) %>%
    group_by(cell,filename) %>%
    summarise(ana_start = quantile(t_ana,probs=t_ana_quantile,na.rm=TRUE)) %>% #find time of first separating sister for start of anaphase
    ungroup() %>%
    mutate(cell=as.integer(cell),
           kittracking_file_str = filename %>% stringr::str_split('/') %>% purrr::map(last)) %>%
    tidyr::unnest()
  
  greta <- inner_join(Data, start_of_anaphase_df,by=c("kittracking_file_str")) %>%
    dplyr::select(-cell.x,-cell.y,-filename.x,-filename.y) %>%
    filter(Time<ana_start-ana_margin)
  
  #for KK distance
  intersister_dist_3d_df <- greta %>%
    group_by(kittracking_file_str,SisterPairID,Time) %>%
    summarise(intersister_dist_3d = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2)) %>%
    #  group_by(kittracking_file_str,SisterPairID) %>%
    #  summarise(av_intersister_dist_3d = median(intersister_dist_3d,na.rm=TRUE)) %>%
    mutate(filename=kittracking_file_str)
  
  thunberg <- greta %>% 
    group_by(kittracking_file_str,SisterPairID,Time) %>%
    summarise(centre_normal_position = mean(Position_1)) %>%
    group_by(kittracking_file_str,SisterPairID) %>%
    mutate(dx=c(NA,diff(centre_normal_position))/dt) %>%
    summarise(centre_normal_speed = sd(dx,na.rm=TRUE)) %>%
    mutate(filename=kittracking_file_str)
  
  periods_df <- greta %>% 
    group_by(kittracking_file_str,SisterPairID, SisterID)  %>%
    mutate(amplitude = 0.5*(caTools::runmax(Position_1,window_size) - 
                              caTools::runmin(Position_1,window_size))) %>%
    summarise(#period = get_period(Position_1,dt),
              amplitude = median(amplitude,na.rm=TRUE)
    )
  long_periods_df <- periods_df %>% 
    tidyr::gather(key=measurement,value=val,-kittracking_file_str,-SisterPairID,-SisterID) %>%
    mutate(filename=kittracking_file_str)
  long_periods_df$measurement <- factor(long_periods_df$measurement,labels = c("Amplitude (um)"))#,"Period (s)"))
  
  
  summarised_draws_df <- draws %>%
    filter(converged) %>%
    group_by(SisterPairID,filename,param) %>%
    summarise(#av_radius = first(av_radius), #calculate av_radius separately so that can include unpaired sisters and short tracks
      theta = median(theta),
      p_reversal=median(p_reversal),
      p_revtoana=median(p_revtoana)) %>%
    add_kittracking_column() %>%
    mutate(param = case_when(param=="p_icoh" ~ "p_coh",
                             param=="p_coh" ~ "p_icoh",
                             TRUE ~ param)) %>%
    tidyr::spread(param,theta) %>%
    ungroup() %>%
    dplyr::select(-filename)
  
  ana_times_df <- draws %>% 
    mutate(t_ana = if_else(param=="t_ana",theta,NA_real_)) %>%
    group_by(cell,filename,SisterPairID) %>%
    summarise(t_ana = median(t_ana,na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(cell=as.integer(cell),
           kittracking_file_str = filename %>% stringr::str_split('/') %>% purrr::map(last)) %>%
    tidyr::unnest()
  
  av_spindle_position_df <- Data %>% inner_join(start_of_anaphase_df,                                              by=c("kittracking_file_str")) %>%
    filter(Time < ana_start-ana_margin) %>%
    group_by(SisterPairID,SisterID,kittracking_file_str) %>%
    summarise(av_spindle_position = (-1)^(1+first(SisterID))*median(Position_1,na.rm=TRUE),
              av_radius = median(sqrt(Position_2^2 + Position_3^2),na.rm=T))
  
  twist_df <- Data %>% inner_join(start_of_anaphase_df,
                                  by=c("kittracking_file_str")) %>%
    filter(Time < ana_start-ana_margin) %>%
    group_by(SisterPairID,kittracking_file_str,Frame) %>%
    summarise(cos_phi = compute_angle_phi(Position_1,Position_2,Position_3)) %>%
    group_by(SisterPairID,kittracking_file_str) %>%
    summarise(twist = acos(median(abs(cos_phi),na.rm=TRUE))/pi*180)
  
  
  asynchrony_df <- ana_times_df %>% group_by(kittracking_file_str) %>%
    summarise(asynchrony = mad(t_ana,na.rm=TRUE))
  
  relative_t_ana_df <- ana_times_df %>% left_join(ana_times_df %>% group_by(kittracking_file_str) %>%
                                                    summarise(t_ana_median = median(t_ana,na.rm=TRUE))) %>%
    mutate(relative_t_ana = t_ana-t_ana_median) %>%
    dplyr::select(-t_ana,-filename,-cell)
  
  all_the_interesting_stats_df <- long_periods_df %>%
    filter(measurement %in% c("Amplitude (um)")) %>% #exclude period, redoing this
    tidyr::spread(measurement,val) %>%
    full_join(intersister_dist_3d_df %>%
                group_by(SisterPairID,kittracking_file_str) %>%
                summarise(intersister_dist_3d = median(intersister_dist_3d,
                                                       na.rm=TRUE)),
              by=c("SisterPairID","kittracking_file_str")) %>%
    full_join(thunberg %>% dplyr::select(-filename),
              by=c("SisterPairID","kittracking_file_str")) %>%
    full_join(summarised_draws_df,by=c("kittracking_file_str","SisterPairID")) %>%
    full_join(av_spindle_position_df,by=c("kittracking_file_str","SisterPairID","SisterID")) %>%
    full_join(twist_df,by=c("kittracking_file_str","SisterPairID")) %>%
    full_join(asynchrony_df,by="kittracking_file_str") %>%
    mutate(tension = kappa*(intersister_dist_3d-L)) %>%
    full_join(relative_t_ana_df,by=c("kittracking_file_str","SisterPairID"))
  
  return(all_the_interesting_stats_df)
} 

generate_figures_based_on_states_and_switching <- function(estimate,sigma_sim,jobset_str,
                                                           identifier,dt=2.05){
  library(ggdist)  
  if (is.null(dim(as.matrix(estimate)))){
  warning('MCMC output does not contain any samples for this cell. Stopping without plotting.\n')
    return(0)
  } #check if contains any samples
  if (length(sigma_sim)==0){
    warning('Hidden states output does not contain samples. Stopping without plotting.\n')
    return(0)
  }
  Data <- process_jobset(jobset_str,K=Inf,max_missing=0.25) %>%
    filter(!is.na(SisterID))
  has_separated_df <- Data %>%
    group_by(SisterPairID,Frame) %>%
    arrange(SisterID,SisterPairID,Frame) %>%
    summarise(kk_dist=-diff(Position_1)) %>%
    group_by(SisterPairID) %>%
    summarise(kk_dist_last = last(kk_dist[!is.na(kk_dist)]),
              separated = kk_dist_last>2.0) #cut off of 2um
  pairIDs <- has_separated_df %>% filter(separated) %>% pull(SisterPairID)
  if (length(pairIDs) != nrow(sigma_sim[[1]])) { pairIDs <- unique(Data$SisterPairID)} #revert to using all pairs if fitted previously with everything
  K <- max(Data$Frame)

  trail_switch_matrix <- purrr::map(sigma_sim, function(x) apply(x,1,count_trail_switches)) %>% do.call(rbind,.)
  lead_switch_matrix <- purrr::map(sigma_sim, function(x) apply(x,1,count_lead_switches)) %>% do.call(rbind,.)
  joint_switch_matrix <- purrr::map(sigma_sim, function(x) apply(x,1,count_joint_switches)) %>% do.call(rbind,.)
  total_switch_matrix <- trail_switch_matrix + lead_switch_matrix + joint_switch_matrix

  proportion_LIDS <- lead_switch_matrix/total_switch_matrix
  proportion_TIDS <- trail_switch_matrix/total_switch_matrix
  
  switching_df_for_cell <- tibble(SisterPairID=pairIDs,
                                  trail=apply(trail_switch_matrix,2,mean),
                                  lead=apply(lead_switch_matrix,2,mean),
                                  joint=apply(joint_switch_matrix,2,mean),
                                  total=apply(total_switch_matrix,2,mean),
                                  proportion_LIDS=apply(proportion_LIDS,2,function(x) mean(x,na.rm=T)),
                                  proportion_TIDS=apply(proportion_TIDS,2,function(x) mean(x,na.rm=T))
                                  )
  LIDS_TIDS_ratio_plt <- switching_df_for_cell %>% 
    rename(Lead=proportion_LIDS,Trail=proportion_TIDS) %>%
    tidyr::gather("switch_type","proportion",Lead,Trail) %>% #proportion_LIDS,proportion_TIDS) %>% 
    ggplot(aes(y=switch_type,x=proportion)) + 
    # stat_halfeye(aes(fill = stat(x < 0.5))) +
    stat_dotsinterval(aes(fill=switch_type),slab_color=NA) +
    theme_bw() +
    theme(legend.position="none") + 
    scale_fill_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4,1)]) + 
    labs(x="Proportion of events",y="Kinetochore initiating switch")
  ggsave(here::here(paste0("plots/LIDS_TIDS_ratio_from_single_cell_",
                           identifier,".eps")),width=4,height=4)

##############
  converged_df <- purrr::map(seq_along(pairIDs), 
                             function(x) check_max_Rhat_less_than_tol_hierarchical(estimate,tol=1.10,trackID=x)) %>%
    bind_rows(.id="track") %>%
    mutate(SisterPairID = pairIDs[as.integer(track)])
if (sum(converged_df$converged)<10){
warning('Fewer than 10 pairs converged in this cell. Estimates not relaible, so stopping early without plotting.\n')
return(0)
}
#results object goes here and contains median parameter estimates for all the key params
#extract parameters for each sister Pair
  results <- tryCatch({
    spread_draws(model=estimate,theta[track,param],p_reversal,p_revtoana) %>% #this is the memory intensive part: xi[track,timept,state]
      summarise(theta = median(theta,na.rm=TRUE),
                xi    = NA, #median(xi,na.rm=TRUE),
                #v_ana = median(v_ana,na.rm=TRUE),
                #t_ana = median(t_ana,na.rm=TRUE),
                p_reversal = median(p_reversal,na.rm=TRUE),
                p_revtoana = median(p_revtoana,na.rm=TRUE))
  }, error = function(err) {
    print(err) #this catches trajectories within a cell when there has been an error in a majority but don't lose the cell
    tibble(param=NA,theta=NA,timept=NA,state=NA,xi=NA,p_reversal=NA,p_revtoana=NA)
  })  %>%
    mutate(SisterPairID = pairIDs[as.integer(track)]) %>%
    left_join(converged_df %>% select(-track),by=c("SisterPairID")) #add info about whether each track converged

jobset_str_list <- jobset_str
draws <- results %>% bind_rows(.id="cell") %>%
  mutate(filename=jobset_str_list[as.integer(cell)]) %>%
  group_by(filename,cell,SisterPairID) %>%
  mutate(param = case_when(
    param==1 ~ "tau",
    param==2 ~ "alpha",
    param==3 ~ "kappa",
    param==4 ~ "v_minus",
    param==5 ~ "v_plus",
    param==6 ~ "p_icoh",
    param==7 ~ "p_coh",
    param==8 ~ "L",
    param==9 ~ "v_ana",
    param==10 ~ "t_ana"
  ))

all_the_interesting_stats_df <- Data %>% 
  mutate(cell=1,filename=jobset_str) %>%
  add_kittracking_column() %>%
  get_all_interesting_stats(draws)

#############

#NB for more than one cell will need to edit this, join is assuming one cell only
all_stats_plus_switching_df <- all_the_interesting_stats_df %>% 
  left_join(switching_df_for_cell,by="SisterPairID")

p1 <- ggplot(all_stats_plus_switching_df,
             aes(av_radius,relative_t_ana,
                 color=lead,size=-v_minus)) + 
  geom_point() + 
  theme_bw() +
  labs(x="Radius within\nmetaphase plate (um)",
       y="Relative time\nof anaphase onset (s)",
       color="# of LIDS") +
  scale_color_gradient(low="#f7fcb9",high="#31a354")
p2 <- ggplot(all_stats_plus_switching_df,
             aes(av_radius,`Amplitude (um)`,
                 color=lead,size=-v_minus)) + 
  geom_point() + 
  theme_bw() +
  labs(x="Radius within\nmetaphase plate (um)",
       y="Amplitude (um)",
       color="# of LIDS") + 
  scale_color_gradient(low="#f7fcb9",high="#31a354")
LIDS_with_multiple_variables_plt <- p1/p2 +  plot_layout(guides = 'collect')
LIDS_with_multiple_variables_plt
ggsave(here::here(paste0("plots/switches_by_plate_position_and_other_variables_",
                         identifier,".eps")),width=4,height=4)

df <- tidyr::gather(all_stats_plus_switching_df,"switch_type","num",trail,lead)

p1 <- ggplot(df,aes(av_radius,num,
                    color=factor(switch_type))) +  
  geom_point() + 
  theme_bw() +
  labs(x="Radius within\nmetaphase plate (um)",
       y="# of switches",color="Switch type") + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4,1)])
p2 <- ggplot(df,aes(`Amplitude (um)`,num,color=factor(switch_type))) + 
  geom_point() + 
  theme_bw() +
  labs(x="Amplitude (um)",y="# of switches",color="Switch type") + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4,1)])
p3 <- ggplot(df,
             aes(relative_t_ana,num,color=factor(switch_type))) + 
  geom_point() + 
  theme_bw() +
  labs(x="Relative time\nof anaphase onset (s)",
       y="# of switches",color="Switch type") + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4,1)])
p4 <- ggplot(df,
             aes(tension,num,color=factor(switch_type))) + 
  geom_point() + 
  theme_bw() +
  labs(x="Tension",y="# of switches",color="Switch type") + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4,1)]) 
p5 <- ggplot(df,
             aes(v_minus,num,color=factor(switch_type))) + 
  geom_point() + 
  theme_bw() +
  labs(x=expression(v["-"] ~ "(um/s)"),y="# of switches",color="Switch type") + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4,1)]) +
  scale_x_continuous(breaks=c(-0.04,-0.06,-0.08))
p6 <- ggplot(df,
             aes(tau,num,color=factor(switch_type))) + 
  geom_point() + 
  theme_bw() +
  labs(x=expression(tau),y="# of switches",color="Switch type") + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(4,1)]) 
switching_and_covariates_plt <- p1+p2+p3+p4+p5+p6+plot_layout(guides = "collect",nrow=3) & theme(legend.position = 'top')
switching_and_covariates_plt
ggsave(here::here(paste0("plots/switches_varying_with_other_variables_two_colour_",
                         identifier,".eps")),width=6,height=4)

get_proportion_of_pairs_in_each_state <- function(mcmc_iter,sigma_sim,K){
  df <- tibble(
    proportion_pp=colMeans(sigma_sim[[mcmc_iter]]==1,na.rm=T),
    proportion_pm=colMeans(sigma_sim[[mcmc_iter]]==2,na.rm=T),
    proportion_mp=colMeans(sigma_sim[[mcmc_iter]]==3,na.rm=T),
    proportion_mm=colMeans(sigma_sim[[mcmc_iter]]==4,na.rm=T),
    proportion_ana=colMeans(sigma_sim[[mcmc_iter]]==5,na.rm=T),
    proportion_rev=colMeans(sigma_sim[[mcmc_iter]]==6,na.rm=T),
    frame=seq_len(K)
  )
}
df <- purrr::map_df(1:length(sigma_sim),
                    function(iter) get_proportion_of_pairs_in_each_state(iter,sigma_sim,K),.id=".iter")

prop_kts_in_each_state_plt <- df %>% 
  tidyr::gather(state,proportion,-frame,-.iter) %>%
  group_by(frame,state) %>%
  summarise(proportion=mean(proportion,na.rm=T)) %>%
  ggplot(aes((frame-1)*dt,proportion,color=state)) + 
  geom_step() + 
  theme_bw() +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(n=4,"PiYG"),
                              "black","grey"),
                     labels=c("A","--","-+","+-","++","R")) +
  labs(x="Time (s)",y="Proportion of KTs")
prop_kts_in_each_state_plt
ggsave(here::here(paste0("plots/proportion_of_kts_in_each_state_",
                         identifier,".eps")),width=4,height=3)

top_left_panel <- LIDS_TIDS_ratio_plt / LIDS_with_multiple_variables_plt + plot_layout(heights=c(1,2))
{top_left_panel | switching_and_covariates_plt} / prop_kts_in_each_state_plt + plot_layout(heights=c(2,1))
ggsave(here::here(paste0("plots/comp_biol_fig6_",identifier,".eps")),
                  width=210,height=297,units = "mm")

av_num_reversal_frames <- purrr::map(seq_along(sigma_sim),function(x) rowSums(sigma_sim[[x]]==6,na.rm=T)) %>% 
  unlist() %>% matrix(nrow=length(sigma_sim),byrow = TRUE) %>% colMeans()
reversals_df <- tibble(SisterPairID=pairIDs,
                       av_num_reversal_frames=av_num_reversal_frames)

num_pairs_with_reversals <- reversals_df %>%
  mutate(has_reversals=av_num_reversal_frames>1) %>%
  group_by(has_reversals) %>% #when more cells with sampled states, group by cell as well
  count() %>%
  tidyr::spread(key=has_reversals,value=n)
if (ncol(num_pairs_with_reversals)==2){
num_pairs_with_reversals_plt <- num_pairs_with_reversals %>%
  ggplot(aes(`TRUE`,`FALSE`)) + 
  geom_point() + 
  geom_abline(slope=-1,intercept=46,linetype="dashed",color="grey") +
  theme_bw() +
  scale_x_continuous(limits=c(0,46)) + 
  scale_y_continuous(limits=c(0,46)) +
  labs(x="# of pairs with reversals",
       y="# of pairs with no reversals")


reversals_vs_sumstats_plt <- all_the_interesting_stats_df %>% 
  inner_join(reversals_df) %>% 
  tidyr::gather(sumstat,value,av_radius,
                v_minus,v_plus,intersister_dist_3d,`Amplitude (um)`,
                tau,relative_t_ana,tension,av_spindle_position) %>% 
  ggplot(aes(value,av_num_reversal_frames)) + 
  geom_point(aes(color=(av_num_reversal_frames>0))) + 
  facet_wrap(.~sumstat,scales="free") +
  theme_bw() +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")[c(1,4)]) +
  labs(x="Summary statistic",y="Expected number of frames in the reversal state") + theme(legend.position="none")
ggsave(here::here(paste0("plots/reversals_vs_other_sumstats_",identifier,".eps")),
       width=210,height=150,units="mm")

reversals_vs_sumstats_plt / {num_pairs_with_reversals_plt | plot_spacer() | plot_spacer() } + plot_layout(heights=c(3,1))
ggsave(here::here(paste0("plots/reversals_fig_",identifier,".eps")),
       width=210,height=150,units="mm")
} else {
cat("skipping reversals plots - check reversals for this cell")
}

return(switching_df_for_cell %>% mutate(filename=jobset_str)) #success
}

