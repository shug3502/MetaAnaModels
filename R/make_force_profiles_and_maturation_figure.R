make_force_profile_and_maturation_figure(){
Data_2s <- purrr::map(jobset_str_list,
                      function(x) process_jobset(x,max_missing=0.25,K=Inf,
                                                 plot_opt=FALSE)) %>%
  bind_rows(.id="cell") %>%
  mutate(filename=jobset_str_list[as.integer(cell)])

draws <- bind_rows(readRDS('~/Documents/Postdoc/Modelling/AnaStanRefactor/AnaStan/fits/median_anaphase_reversals_parameter_estimates_JHprocess_LIDS_TIDS_2s_v303.rds'),
                   readRDS('~/Documents/Postdoc/Modelling/AnaStanRefactor/AnaStan/fits/median_anaphase_reversals_parameter_estimates_JHprocess_LIDS_TIDS_2s_v304.rds'))

ktf <- Data_2s %>% add_kittracking_column() %>% pull(kittracking_file_str) %>% first()
Data <- Data_2s %>% add_kittracking_column() %>% filter(kittracking_file_str==ktf)
pairIDs <- Data %>% pull(SisterPairID) %>% unique()
K=307
get_states_df <- function(mcmc_iter,sigma_sim,pairIDs,K){
  states_df <- tibble(state=c(t(sigma_sim[[mcmc_iter]])),
                      SisterPairID=rep(pairIDs,K),
                      Frame=rep(seq_len(K),length(pairIDs)))
}
states_df <- purrr::map_df(1:4,function(x) get_states_df(x,sigma_sim,pairIDs,K),.id="iter")

forces_df <- states_df %>% mutate(kittracking_file_str=ktf) %>%
  full_join(ana_times_df,by=c("kittracking_file_str","SisterPairID")) %>%
  mutate(Frames_since_start = Frame - t_ana%/%dt) %>%
  dplyr::select(-t_ana) %>%
  inner_join(complex_positions_by_frame_treatment_df,by=c("kittracking_file_str","SisterPairID","Frames_since_start")) %>%
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
}