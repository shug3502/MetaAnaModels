make_anaphase_times_and_speed_figure <- function(jobset_str_list,draws){

Data_2s <- purrr::map(jobset_str_list,
                      function(x) process_jobset(x,max_missing=0.25,K=Inf,
                                                 plot_opt=FALSE)) %>%
  bind_rows(.id="cell") %>%
  mutate(filename=jobset_str_list[as.integer(cell)]) %>%
  add_kittracking_column()

  ana_margin <- 60
  ana_times_df <- draws %>% 
    mutate(t_ana = if_else(param=="t_ana",theta,NA_real_)) %>%
    group_by(cell,filename,SisterPairID) %>%
    summarise(t_ana = median(t_ana,na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(cell=as.integer(cell),
           kittracking_file_str = filename %>% stringr::str_split('/') %>% purrr::map(last)) %>%
    tidyr::unnest()
  
  start_of_anaphase_df <- draws %>%
    mutate(t_ana = if_else(param=="t_ana",theta,NA_real_)) %>%
    group_by(cell,filename) %>%
    summarise(ana_start = quantile(t_ana,probs=0.5,na.rm=TRUE)) %>% #find time of first separating sister for start of anaphase
    ungroup() %>%
    mutate(cell=as.integer(cell),
           kittracking_file_str = filename %>% stringr::str_split('/') %>% purrr::map(last)) %>%
    tidyr::unnest()
  
  av_spindle_position_df <- Data_2s %>% inner_join(start_of_anaphase_df,                                              
                                                   by=c("kittracking_file_str")) %>%
    filter(Time < ana_start-ana_margin) %>%
    group_by(SisterPairID,SisterID,kittracking_file_str) %>%
    summarise(av_spindle_position = (-1)^(1+first(SisterID))*median(Position_1,na.rm=TRUE),
              av_radius = median(sqrt(Position_2^2 + Position_3^2),na.rm=T))
  
  relative_t_ana_df <- ana_times_df %>% left_join(ana_times_df %>% group_by(kittracking_file_str) %>%
                                                    summarise(t_ana_median = median(t_ana,na.rm=TRUE))) %>%
    mutate(relative_t_ana = t_ana-t_ana_median) %>%
    dplyr::select(-t_ana,-filename,-cell)
  
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
  
  simple_positions_by_frame_treatment_df <- Data_2s %>%
    full_join(start_of_anaphase_df,by="kittracking_file_str") %>%
    mutate(Frames_since_start = Frame - ana_start%/%dt) %>%
    filter(Frames_since_start < nframes_early) %>%
    filter(Frames_since_start >= -nframes_early) %>%
    group_by(kittracking_file_str,Frames_since_start,SisterPairID,SisterID) %>%
    summarise(Position_1 = first(Position_1),
              Position_2 = first(Position_2),
              Position_3 = first(Position_3)) %>%
    group_by(kittracking_file_str) %>%
    filter(length(unique(SisterPairID))>min_num_sisters)
  
green_val <- RColorBrewer::brewer.pal(4,"PiYG")[4]
t_ana_df <- inner_join(av_spindle_position_df,relative_t_ana_df,by=c("kittracking_file_str","SisterPairID"))
mdl <- lm(relative_t_ana ~ av_radius, data=t_ana_df)
summary(mdl) %>% print()

plt1 <- ggplot(t_ana_df,
               aes(av_radius,relative_t_ana)) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=FALSE,color=green_val) + 
  theme_bw() + 
  scale_x_continuous(limits=c(0,6)) +
  scale_y_continuous(limits=c(-100,100)) + 
  labs(x="Average radius (um)\nduring metaphase",
       y="Anaphase onset time (s)\nof each pair relative\nto median for the cell",
       title="KT pairs at the edge\n of the metaphase plate\nseparate earlier"
  )

v_ana_1d_df <- inner_join(av_spindle_position_df,
           summarised_draws_df,by=c("kittracking_file_str","SisterPairID"))
mdl <- lm(v_ana ~ av_radius,data=v_ana_1d_df)
summary(mdl) %>% print()
plt2 <- ggplot(v_ana_1d_df, 
               aes(av_radius,v_ana)) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=FALSE,color=green_val) + 
  scale_x_continuous(limits=c(0,6)) +
  theme_bw() + 
  #  scale_y_continuous(limits=c(-100,100)) + 
  labs(x="Average radius (um)\n during metaphase",
       y="Anaphase speed (um/s)",
       title="KT pairs at the edge of\nthe metaphase plate\nhave lower speed v_A"
  )

framewise_speed_df <- simple_positions_by_frame_treatment_df %>% 
  filter(Frames_since_start > 0) %>% 
  group_by(kittracking_file_str,SisterPairID,SisterID) %>%
  mutate(framewise_dist_travelled=c(0,sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2))) %>%
  ungroup() %>% drop_na()%>%
  group_by(kittracking_file_str,SisterPairID,SisterID) %>%
  summarise(speed_3d=median(framewise_dist_travelled/dt,na.rm=T)) %>%
  inner_join(av_spindle_position_df,by=c("kittracking_file_str", "SisterPairID", "SisterID")) %>%
  filter(av_radius<=8)

mdl <- lm(speed_3d ~ av_radius,data=framewise_speed_df)
summary(mdl) %>% print()

plt3 <- ggplot(framewise_speed_df,aes(av_radius,speed_3d)) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=FALSE,color=green_val) + 
  theme_bw() + 
  labs(x="Average radius (um)\nduring metaphase",
       y="Average 3D framewise speed\nin anaphase (um/s)",
       title="3D framewise speed in\nanaphase does not show\nsame spatial effect")

plt1 | plt2 | plt3 
ggsave(here::here("plots/anaphase_time_and_speed_versus_radius.eps"), width=210,height=150,units="mm")
}

