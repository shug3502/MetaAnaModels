make_population_params_and_heterogeneity_figure <- function(){
jobset_str_list <- list.files(path="data/",pattern="*.csv",
                              recursive = TRUE,full.names = TRUE)
Data_2s <- purrr::map(jobset_str_list,
                      function(x) process_jobset(x,max_missing=0.95,K=Inf,
                                                 plot_opt=FALSE)) %>%
  bind_rows(.id="cell") %>%
  mutate(filename=jobset_str_list[as.integer(cell)])

draws <- bind_rows(readRDS('~/Documents/Postdoc/Modelling/AnaStanRefactor/AnaStan/fits/median_anaphase_reversals_parameter_estimates_JHprocess_LIDS_TIDS_2s_v303.rds'),
                   readRDS('~/Documents/Postdoc/Modelling/AnaStanRefactor/AnaStan/fits/median_anaphase_reversals_parameter_estimates_JHprocess_LIDS_TIDS_2s_v304.rds'))

individual_plt <- draws %>%
  filter(param %in% c("tau","alpha","kappa","v_minus","v_plus","L")) %>%
  ggplot(aes(theta)) + 
  geom_histogram(aes(y=..count../sum(..count..))) + 
  facet_wrap(.~param,labeller=label_parsed,scales="free_x") + 
  theme_bw() +
  labs(x="Parameter",y="Proportion of kinetochore pairs")

shared_plt <- draws %>%
  filter(param %in% c("p_icoh","p_coh")) %>%
  tidyr::spread(param,theta) %>%
  rename(p_icoh=p_coh,p_coh=p_icoh,p_AR=p_reversal,p_RA=p_revtoana) %>%
  tidyr::gather("param","theta",p_icoh,p_coh,p_AR,p_RA) %>%
  group_by(filename,cell,param) %>%
  summarise(theta=median(theta)) %>%
  ggplot(aes(theta)) + 
  geom_histogram(aes(y=..count../sum(..count..))) + 
  facet_wrap(.~param,labeller=label_parsed,scales="free_x") + 
  theme_bw() +
  labs(x="Parameter",y="Proportion of cells")

draws_with_positions_df <- Data_2s %>% add_kittracking_column() %>% 
  inner_join(draws %>% add_kittracking_column(),by=c("kittracking_file_str","SisterPairID"))

a <- draws_with_positions_df %>%
  mutate(r = cut(av_radius,c(0,2,4,6))) %>%
  rbind(draws_with_positions_df %>% mutate(r="All")) %>%
  group_by(kittracking_file_str,param,r) %>% summarise(within_cell_spread=sd(theta,na.rm=T)) %>% 
  filter(param %in% c("tau","alpha","kappa","v_minus","v_plus","L")) %>%
  ungroup() %>%
  mutate(r=factor(r,levels=c("All","(0,2]","(2,4]","(4,6]")))
b <- draws_with_positions_df %>%
  mutate(r = cut(av_radius,c(0,2,4,6))) %>%
  rbind(draws_with_positions_df %>% mutate(r="All")) %>%
  group_by(kittracking_file_str,param,r) %>% 
  summarise(theta=median(theta,na.rm=T)) %>%
  group_by(param,r) %>% 
  summarise(between_cell_spread=sd(theta,na.rm=T)) %>%
  filter(param %in% c("tau","alpha","kappa","v_minus","v_plus","L")) %>%
  ungroup() %>%
  mutate(r=factor(r,levels=c("All","(0,2]","(2,4]","(4,6]")))

h1 <- ggplot(a,aes(r,within_cell_spread)) + 
  geom_violin(draw_quantiles = 0.5) + 
  geom_jitter() +
  facet_wrap(.~param,labeller=label_parsed,scales="free") + 
  theme_bw()
h2 <- h1 + 
  geom_point(data=b,
             aes(x=r, y = between_cell_spread),color="red",size=4) + 
  labs(x="Radius in metaphase plate (um)",y="Standard deviation of parameter")

individual_plt / shared_plt / h2 + plot_layout(heights=c(1,1,1))
ggsave(here::here("plots/comp_biol_fig5.eps"),device=cairo_ps,width=210,height=297,units="mm")

#Assuming a normal distribution for within cell spread, evaluate probability that between cell spread comes from this distribution
a %>% ungroup() %>% filter(abs(within_cell_spread)>0) %>% 
  group_by(param,r) %>% 
  summarise(mu=mean(within_cell_spread),
            sig=sd(within_cell_spread)) %>% 
  inner_join(b) %>% 
  mutate(pval=pnorm(between_cell_spread,mean=mu,sd=sig,lower.tail=TRUE))
}