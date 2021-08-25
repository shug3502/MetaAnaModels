args = commandArgs(trailingOnly=TRUE)
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

rstan::rstan_options(auto_write = TRUE) #tries to avoid recompiling stan code
options(mc.cores = parallel::detectCores()) #uses as many cores as you have
source(here::here('R/helper_fns.R'))
source(here::here('R/fit_anaphase_changept_model.R'))
source(here::here('R/fit_anaphase_reversals_model.R'))
source(here::here('R/run_anaphase_models_for_ith_jobset.R'))
source(here::here('R/extract_hidden_states.R'))
##############################
#options etc
dt=2.05
run_analysis <- TRUE
use_parallel <- FALSE
run_changept_anyway <- FALSE
num_iter <- 100
identifier = args[1]
path_to_folder = args[2]

print(identifier)
print(path_to_folder)
#############################

jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",full.names=TRUE)
K_list <- rep(NA,length(jobset_str_list))
stopifnot(length(jobset_str_list)>0)

#jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame//kittracking001-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture10_flowdec_deconvolved.ome.csv")
#data_single_pair <- process_jobset(jobset_str,K=Inf,max_missing=0.25) %>%
#  filter(!is.na(SisterID)) #omit unpaired KTs
#g <- ggplot(data_single_pair, aes(x=Time, y=Position_1,color=factor(SisterID))) +
#  geom_line() +
#  facet_wrap(.~SisterPairID) +
#  theme_bw() +
#  theme(legend.position = "None",
#        strip.text.x = element_text(size = 8),
#        axis.text.x = element_text(angle=90)) +
#  labs(x="Time (s)",y="Position (um)")
#g
#ggsave(here::here("plots/200818_untreated_capture10_all_pairs.eps"),width=7,height=7)

for (i in seq_along(jobset_str_list)){
  #fit changepoint model and hierarchical model
  out <- run_anaphase_models_for_ith_jobset(i,jobset_str_list[i],K_list[i],here::here("fits"),identifier,run_analysis,use_parallel,run_changept_anyway,dt,num_iter)
}  
#extract hidden states
#sigma_sim <- extract_hidden_states(out)

