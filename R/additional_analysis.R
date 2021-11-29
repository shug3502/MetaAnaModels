library(rstan) 
library(dplyr) 
library(purrr)
library(zoo) 
library(here) 
library(ggplot2) 
library(patchwork) 
library(bayesplot)
library(tidybayes)
library(stringr)
library(readr)
library(tidyr)
source('R/extract_hidden_states.R')
source('R/helper_fns.R')
#options etc
dt=2.05
load_draws <- TRUE
fits_folder_str <- "fits"
identifier = "jonathanharrison_MetaAnaModels_2s_v411" #args[1]
path_to_folder = "data" #args[2]

jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",
                              full.names=TRUE,recursive=TRUE)
stopifnot(length(jobset_str_list)>0)
#jobset_str_list <- jobset_str_list[1:15]
if (load_draws){
  draws <- readRDS('fits/median_anaphase_reversals_parameter_estimates_jonathanharrison_MetaAnaModels_2s_v411.rds')
} else {
  draws <- summarise_batch_anaphase_reversals(identifier,jobset_str_list,dt=dt,
                                     fits_folder_str = fits_folder_str,
                                     max_missing = 0.25,tol = 1.05,
                                     run_analysis = FALSE,
                                     use_parallel = FALSE)
}
jobset_str_list <- draws %>% filter(converged) %>% pull(filename) %>% unique() #only use cells where MCMC has converged
draws <- draws %>% filter(converged)

source('R/how_common_are_reversals.R')
reversals_per_pair <- tibble()
radius_per_pair <- tibble()
laziness_per_pair <- tibble()
for (i in seq_along(jobset_str_list)){
  #extract hidden states
#  jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame//kittracking001-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture10_flowdec_deconvolved.ome.csv")
  job_id = stringr::str_split(jobset_str_list[i],"kittracking")[[1]][2]
  edited_job_id = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  path_to_est <- file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',
                                                 edited_job_id,'.rds',sep=''))
  if (file.exists(path_to_est)){
    estimate <- readRDS(file=path_to_est)
    sigma_sim <- extract_hidden_states(estimate)
    rr <- how_common_are_reversals(jobset_str_list[i],sigma_sim)
    reversals_per_pair <- bind_rows(reversals_per_pair,rr)
    radius_per_pair <- bind_rows(radius_per_pair,get_av_radius_of_pairs(jobset_str_list[i]))
    laziness_per_pair <- bind_rows(laziness_per_pair,get_laziness_for_cell(jobset_str_list[i],draws,min_num_sisters=10))    
    success=1
  } else {
    cat('Files from sampling do not exist. Skipping this cell.\n')
    success=0
  }
  if (success){cat("SUCCESS!\n")} else {cat("FAILED :(\n")}
}

g <- ggplot(reversals_per_pair, aes(reversals)) + geom_histogram() + theme_bw() + labs(x="Average number of frames in the reversal state",y="Number of kinetochore pairs")
#g
#ggsave(here::here("plots/histogram_of_frames_in_reversal_state.eps"),width=6,height=4)

reversals_per_pair %>% mutate(long_reversal=reversals>10) %>% count(long_reversal)

h <- ggplot(reversals_per_pair %>% mutate(long_reversal=reversals>10) %>% 
  group_by(filename) %>%
  count(long_reversal) %>% 
  group_by(filename) %>%
  summarise(total=sum(n),
            proportion=1-first(n)/total),
  aes(proportion)) + 
  geom_histogram() + 
  theme_bw() +
  labs(x="Proportion of kinetochore pairs\nwith >10 frames in reversal state",
       y="Number of cells") 
#h
#ggsave(here::here("plots/histogram_of_cells_with_long_reversals.eps"),width=6,height=4)

g | h
ggsave(here::here("plots/histogram_reversals_both.eps"),width=210,height=100,units="mm")

lots_of_info_df <- draws %>%
inner_join(reversals_per_pair,by=c("filename","SisterPairID")) %>%
#inner_join(radius_per_pair,by=c("filename","SisterPairID")) %>%
inner_join(laziness_per_pair,by=c("filename","SisterPairID"))

p1 <- lots_of_info_df %>% 
  mutate(kk_dist_discrete=cut(av_intersister_dist,breaks=seq(from=0.8,to=1.6,by=0.2))) %>%
  ggplot(aes(kk_dist_discrete,reversals)) +
    geom_jitter() + 
    geom_violin(draw_quantiles=0.5) +
    theme_bw() +
    labs(x="Intersister distance (um)",y="Average number of frames in the reversal state") 
p1 
ggsave(here::here("plots/kk_dist_reversals_violin.eps"),width=210,height=100,units="mm")  

p2 <- lots_of_info_df %>% 
  mutate(reversals_discrete=cut(reversals,breaks=c(0,10,30))) %>%
  ggplot(aes(reversals_discrete,max_laziness)) +
    geom_jitter() +
    geom_violin(draw_quantiles=0.5) +
    theme_bw() +
    labs(x="Laziness",y="Average number of frames in the reversal state")
p2
ggsave(here::here("plots/laziness_reversals_violin.eps"),width=210,height=100,units="mm")

library(ggpubr)
library(rstatix)
stat.test <- lots_of_info_df %>%
      mutate(laziness_discrete=cut(max_laziness,breaks=c(-Inf,1.93,3))) %>%
      filter(!is.na(theta)) %>%
      group_by(param) %>%
      wilcox_test(theta ~ laziness_discrete) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    stat.test2 <- stat.test %>%
      add_xy_position(step.increase=0) %>%
      group_by(interaction(group1,group2)) %>%
      mutate(ind=group_indices(),
             y.position=y.position*(1+0.06*ind),
             p.adj=p_format(p.adj,accuracy=10^-3,digits=1,add.p=TRUE)) %>%
      ungroup()

p3 <- lots_of_info_df %>%
  filter(param!="t_ana") %>%
  mutate(laziness_discrete=cut(max_laziness,breaks=c(-Inf,2,3))) %>%
  ggplot(aes(laziness_discrete,theta)) +
    geom_jitter(aes(color=laziness_discrete)) +
    geom_violin(draw_quantiles=0.5) +
    scale_colour_manual(values = RColorBrewer::brewer.pal(4,"PiYG")[c(1,4)]) +
    facet_wrap(.~param,scales="free",ncol=3,labeller=label_parsed) +
    theme_bw() +
    labs(x="Max laziness",
         y="Parameter") +
    theme(legend.position = "None") +
#          axis.text.x = element_text(angle = -45,vjust=0)) +
  stat_pvalue_manual(stat.test2,
                     hide.ns = TRUE,label = "{p.adj}",size=2
) +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.10)))
p3
ggsave(here::here("plots/laziness_reversals_violin.eps"),width=210,height=100,units="mm")






