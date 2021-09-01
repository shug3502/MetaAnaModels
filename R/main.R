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
##############################
#options etc
dt=2.05
run_analysis <- TRUE
use_parallel <- FALSE
run_changept_anyway <- FALSE
num_iter <- 500
identifier = args[1]
path_to_folder = args[2]
fits_folder_str <- here::here("fits")

print(identifier)
print(path_to_folder)
#############################

jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",full.names=TRUE)
K_list <- rep(Inf,length(jobset_str_list))
nJobsets <- length(jobset_str_list)
stopifnot(nJobsets>0)

#for (i in seq_along(jobset_str_list)){
for (i in sample.int(nJobsets,size=nJobsets,replace=FALSE)){ #go through jobsets in random orde
  job_id = stringr::str_split(jobset_str_list[i],"kittracking")[[1]][2]
  output_exists <- (file.exists(file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',
                  paste(job_id %>% stringr::str_replace_all("\\.",""),identifier,sep=""),
                  '.rds',sep='')))) #no need to run again if output already exists
  run <- (run_analysis & !output_exists)
  #fit changepoint model and hierarchical model
  out <- run_anaphase_models_for_ith_jobset(i,jobset_str_list,K_list,fits_folder_str,identifier,run,use_parallel,run_changept_anyway,dt,num_iter)
}  

