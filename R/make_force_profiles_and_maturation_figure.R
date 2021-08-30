make_force_profile_and_maturation_figure <- function(jobset_str_list,identifier,
                                         num_iter=200,dt=2.05,min_num_sisters=10,nframes_early=300,
                                         fits_folder_str="fits"){
  #warning on memory usage for all MCMC iterations across many cells
Data_2s <- purrr::map(jobset_str_list,
                      function(x) process_jobset(x,max_missing=0.25,K=Inf,
                                                 plot_opt=FALSE)) %>%
  bind_rows(.id="cell") %>%
  mutate(filename=jobset_str_list[as.integer(cell)]) %>%
  add_kittracking_column()

draws <- bind_rows(readRDS('~/Documents/Postdoc/Modelling/AnaStanRefactor/AnaStan/fits/median_anaphase_reversals_parameter_estimates_JHprocess_LIDS_TIDS_2s_v303.rds'),
                   readRDS('~/Documents/Postdoc/Modelling/AnaStanRefactor/AnaStan/fits/median_anaphase_reversals_parameter_estimates_JHprocess_LIDS_TIDS_2s_v304.rds'))
ktf_all <- unique(Data_2s$kittracking_file_str) #ids of all cells

#generate ana times df and complex positions etc
ana_times_df <- draws %>%
  mutate(t_ana = if_else(param=="t_ana",theta,NA_real_)) %>%
  group_by(cell,filename,SisterPairID) %>%
  summarise(t_ana = median(t_ana,na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(cell=as.integer(cell),
         kittracking_file_str = filename %>% stringr::str_split('/') %>% purrr::map(last)) %>%
  tidyr::unnest()

#align to anaphase for each trajectory rather than the median for each cell
complex_positions_by_frame_treatment_df <- Data_2s %>%
  full_join(ana_times_df,by=c("kittracking_file_str","SisterPairID")) %>%
  mutate(Frames_since_start = Frame - t_ana%/%dt) %>%
  filter(Frames_since_start < nframes_early) %>%
  filter(Frames_since_start >= -nframes_early) %>%
  group_by(kittracking_file_str,Frames_since_start,SisterPairID,SisterID) %>%
  summarise(Position_1 = first(Position_1),
            Position_2 = first(Position_2),
            Position_3 = first(Position_3)) %>%
  group_by(kittracking_file_str) %>%
  filter(length(unique(SisterPairID))>min_num_sisters) %>%
  mutate(filename=kittracking_file_str)

get_states_df <- function(mcmc_iter,sigma_sim,pairIDs,K){
  states_df <- tibble(state=c(t(sigma_sim[[mcmc_iter]])),
                      SisterPairID=purrr::map(pairIDs,function(x) rep(x,K)) %>% unlist(),
                      Frame=rep(seq_len(K),length(pairIDs)))
}
get_states_df_for_cell <- function(ktf,Data_2s,num_iter){
  Data <- Data_2s %>% filter(kittracking_file_str==ktf)
  jobset_str <- Data$filename
  pairIDs <- Data %>% pull(SisterPairID) %>% unique()
  K <- max(Data$Frame)
  sigma_sim <- extract_hidden_states(NULL,jobset_str,fits_folder_str=fits_folder_str,identifier=identifier) #load sigma_sim
  states_df <- purrr::map_df(1:num_iter,function(x) get_states_df(x,sigma_sim,pairIDs,K),
                           .id="iter")
}

all_states_df <- purrr::map_df(ktf_all,function(x) get_states_df_for_cell(x,Data_2s,num_iter),
                               .id="cell") %>%
  mutate(kittracking_file_str=ktf_all[as.integer(cell)])
  
forces_df <- all_states_df %>%
  inner_join(ana_times_df,by=c("kittracking_file_str","SisterPairID")) %>%
  mutate(Frames_since_start = Frame - t_ana%/%dt) %>%
  dplyr::select(-t_ana,-cell.x,-cell.y,-Frame,-filename) %>%
  inner_join(complex_positions_by_frame_treatment_df,
             by=c("kittracking_file_str","SisterPairID","Frames_since_start")) %>%
  inner_join(draws %>% dplyr::select(-state,-timept) %>% 
               tidyr::spread(key=param,value=theta) %>% 
               add_kittracking_column(),by=c("kittracking_file_str","SisterPairID")) %>% 
  group_by(iter,Frames_since_start,kittracking_file_str,SisterPairID) %>%
  summarise(spring = first(kappa)*(-diff(Position_1)-first(L))*(state<5),
            PEF = first(alpha)*max(Position_1)*(state<5),
            `K-fibre (Meta)` = case_when(state==1 ~ first(v_plus),
                                         state==2 ~ first(v_plus),
                                         state==3 ~ -first(v_minus),
                                         state==4 ~ -first(v_minus),
                                         state==5 ~ 0,
                                         state==6 ~ 0,
                                         TRUE ~ NA_real_),
            `K-fibre (Ana)` = case_when(state==1 ~ 0,
                                        state==2 ~ 0,
                                        state==3 ~ 0,
                                        state==4 ~ 0,
                                        state==5 ~ first(v_ana),
                                        state==6 ~ 0,
                                        TRUE ~ NA_real_))

plt1 <- forces_df %>%
  tidyr::gather(key=force,value=value,-kittracking_file_str,-SisterPairID,-Frames_since_start, -iter) %>%
  group_by(Frames_since_start,kittracking_file_str,SisterPairID,force) %>%
  summarise(value=mean(value,na.rm=T)) %>%
  group_by(force,Frames_since_start) %>%
  summarise(med=median(value,na.rm=T),
            ub=quantile(value,0.9,na.rm=T),
            lb=quantile(value,0.1,na.rm=T)) %>%
  ungroup() %>%
  filter(Frames_since_start<30) %>% filter(Frames_since_start > -60) %>%
  ggplot(aes(x=Frames_since_start*dt,y=med,color=factor(force))) + 
  geom_line() +
  scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")) +
  theme_bw() +
  theme(legend.position = "top") +
  guides(color=guide_legend(ncol=2)) +
  labs(x="Time (s)",y="Force (um/s)",color="")
plt1
ggsave(here::here("plots/force_profiles_around_anaphase.eps"),device=cairo_ps,width=4.5,height=4.5)

plt2 <- forces_df %>% mutate(is_ana=factor(if_else(`K-fibre (Ana)`>0,"Ana","Meta"),levels=c("Meta","Ana"))) %>% 
  tidyr::gather(key=force,value=value,
                -kittracking_file_str,-SisterPairID,
                -Frames_since_start, -iter, -is_ana) %>% 
  group_by(force,is_ana) %>% 
  summarise(val = mean(value,na.rm=T)) %>% 
  drop_na() %>% 
  ggplot(aes(is_ana,val,fill=force)) + 
  geom_bar(stat="identity") + 
  labs(x=" ",y="Force (um/s)" ,fill=" ") + 
  theme_bw() +
  theme(legend.position = "top") +
  guides(fill=guide_legend(ncol=2)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG"))

plt1 | plt2
ggsave(here::here("plots/force_profiles_around_anaphase_and_barchart.eps"),device=cairo_ps,width=4.5,height=4.5)

#######
#maturation in metaphase


simple_positions_by_frame_treatment_df <- Data_2s %>%
  add_kittracking_column() %>%
  full_join(start_of_anaphase_df,by="kittracking_file_str") %>%
  mutate(Frames_since_start = Frame - ana_start%/%dt) %>%
  filter(Frames_since_start < nframes_early) %>%
  filter(Frames_since_start >= -nframes_early) %>%
  group_by(kittracking_file_str,Frames_since_start,SisterPairID,SisterID) %>%
  summarise(Position_1 = first(Position_1),
            Position_2 = first(Position_2),
            Position_3 = first(Position_3)) %>%
  group_by(kittracking_file_str) %>%
  filter(length(unique(SisterPairID))>min_num_sisters) %>%
  mutate(filename=kittracking_file_str) %>% add_treatment_column()

#plate width plot
plate_width_df <- simple_positions_by_frame_treatment_df %>%
  filter(Frames_since_start<=0) %>%
  group_by(kittracking_file_str,Frames_since_start) %>%
  summarise(plate_width=mad(Position_1,na.rm=T)) %>%
  group_by(Frames_since_start) %>%
  summarise(med=median(plate_width,na.rm=T),ub=quantile(plate_width,0.9),lb=quantile(plate_width,0.1))
h1 <- plate_width_df %>%
  ggplot(aes(x=Frames_since_start*dt)) +
  geom_line(aes(y=med)) +
  geom_ribbon(aes(ymin=lb,ymax=ub),color="grey",alpha=0.6) + theme_bw() +
  labs(x="Time (s)",y="Plate width (um)")
h1
# ggsave(here::here("plots/metaphase_plate_width_over_time_2s_data.eps"),device=cairo_ps,width=6,height=4)

twist_df <- complex_positions_by_frame_treatment_df %>%   
  group_by(SisterPairID,kittracking_file_str,Frames_since_start) %>%
  summarise(cos_phi = compute_angle_phi(Position_1,Position_2,Position_3),
            twist = acos(abs(cos_phi))/pi*180)
j1 <- twist_df %>%
  filter(Frames_since_start<=0) %>%
  group_by(Frames_since_start) %>%
  summarise(med=median(twist,na.rm=T),
            ub=quantile(twist,0.9,na.rm=T),
            lb=quantile(twist,0.1,na.rm=T)) %>%
  ggplot(aes(x=Frames_since_start*dt)) + 
  geom_line(aes(y=med)) +
  geom_ribbon(aes(ymin=lb,ymax=ub),color="grey",alpha=0.6) + 
  theme_bw() +
  labs(x="Time (s)",
       y="Twist (degrees)")
# ggsave(here::here("plots/twist_over_time_2s_data.eps"),device=cairo_ps,width=6,height=4)

kk_dist_df <- complex_positions_by_frame_treatment_df %>%   
  group_by(SisterPairID,kittracking_file_str,Frames_since_start) %>%
  summarise(kk_dist = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2))

k1 <- kk_dist_df %>%
  filter(Frames_since_start<=0) %>%
  group_by(Frames_since_start) %>%
  summarise(med=median(kk_dist,na.rm=T),
            ub=quantile(kk_dist,0.9,na.rm=T),
            lb=quantile(kk_dist,0.1,na.rm=T)) %>%
  ggplot(aes(x=Frames_since_start*dt)) + 
  geom_line(aes(y=med)) +
  geom_ribbon(aes(ymin=lb,ymax=ub),color="grey",alpha=0.6) + 
  theme_bw() +
  labs(x="Time (s)",
       y="Intersister distance (um)")
# ggsave(here::here("plots/kk_dist_over_time_2s_data.eps"),device=cairo_ps,width=6,height=4)

h1 | j1 | k1
ggsave(here::here("plots/maturation_over_time_2s_data.eps"),device=cairo_ps,width=6,height=4)

} 
