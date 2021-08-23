args = commandArgs(trailingOnly=TRUE)
devtools::load_all(".")

if (length(args)==0) {
  cat("No inputs supplied. Assuming defaults. Check you intended this\n")
  identifier <- "JHnerov621"
  path_to_folder = "../Nocodazole_washout/"
  dt=4.7
} else if (length(args)>=1) {
  identifier = args[1]
  path_to_folder = args[2]
  dt = as.numeric(args[3])
}
print(identifier)
print(path_to_folder)
print(dt)

#folders <- stringr::str_split(path_to_folder,"/") #convert to path on nero versus path in petabyte
#local_folder <- folders[[1]][length(folders[[1]])]
#path_to_folder <- paste0("../reprocessing_nocodazole_washout/",path_to_folder)
path_to_folder <- paste0("~/data/",path_to_folder)

jobset_str_list <- list.files(path = path_to_folder,pattern="\\.csv$",full.names=TRUE)
cat(paste0("Found ", length(jobset_str_list), " files in the folder ", path_to_folder, "\n"))
fits_folder_str <- here::here("fits")
if (is.na(dt)){ dt=4.7 }
K_list <- rep(Inf,length(jobset_str_list)) #try to use all the movies, but know this is a bad idea in some cases
run_analysis <- TRUE
use_parallel <- FALSE
run_changept_anyway <- FALSE #run_analysis #due to dt problems

#use filename annotations to determine if need to exclude as does not contain anaphase
is_metaphase_only <- purrr::map(jobset_str_list, function(x) stringr::str_detect(x,"_M.ome")) %>% as.logical
jobset_str_list <- jobset_str_list[!is_metaphase_only]
K_list <- K_list[!is_metaphase_only]

run_anaphase_models_for_ith_jobset <- function(i,jobset_str_list,K_list,fits_folder_str,identifier,run_analysis,use_parallel,dt){
  cat(paste("Analysing jobset number ",i,": ",jobset_str_list[i],"\n",sep=""))
  job_id = stringr::str_split(jobset_str_list[i],"kittracking")[[1]][2]
if (file.exists(file.path(fits_folder_str,paste('anaphase_changepoint_estimates',
		paste(job_id %>% stringr::str_replace_all("\\.",""),identifier,sep=""),
		'_change_pt.rds',sep='')))){
cat(file.path(fits_folder_str,paste('anaphase_changepoint_estimates',
                paste(job_id %>% stringr::str_replace_all("\\.",""),identifier,sep=""),
                '_change_pt.rds',sep='')))
cat("(changept) file exists so no need to run again? \n")
run_changept <- run_changept_anyway
} else {
cat("(changept) file does not exists so now attempting to run \n")
run_changept <- TRUE
}
tryCatch({
changept_out <- fit_anaphase_changept_model(jobset_str_list[i], K=K_list[i],
                                            identifier = paste(identifier,"_change_pt",sep=""),
                                            run_analysis = run_changept,
                                            fits_folder_str = fits_folder_str, dt=dt)
t_ana_input_df <- changept_out %>%
  spread_draws(t_ana[track]) %>%
  summarise(t_ana=median(t_ana))

out <- fit_anaphase_reversals_model(jobset_str_list[i],t_ana_input_df, K=K_list[i],
                                           identifier = identifier,
                                           run_analysis = run_analysis,
                                           use_parallel = use_parallel,
                                           fits_folder_str = fits_folder_str,
					   plot_opt=1, dt=dt,
					   num_iter=200 
)
cat("Done ...\n")
},
error = function (err) { print(err)
cat("Continuing but analysis not done for that job ... \n")
})
return(0)
}

if (use_parallel){
library(furrr)
plan(multiprocess,workers=parallel::detectCores())
z_ignore <- furrr::future_map(seq_along(jobset_str_list),function(x) run_anaphase_models_for_ith_jobset(x,
    jobset_str_list,K_list,fits_folder_str,identifier,run_analysis,use_parallel,dt))
} else {
for (ii in seq_along(jobset_str_list)){
  z_ignore <- run_anaphase_models_for_ith_jobset(ii,jobset_str_list,K_list,fits_folder_str,identifier,run_analysis,use_parallel,dt)
}
#z_ignore <- purrr::map(seq_along(jobset_str_list),function(x) run_anaphase_models_for_ith_jobset(x,
#    jobset_str_list,K_list,fits_folder_str,identifier,run_analysis,use_parallel,dt))
}

##use cell newly defined function
cat("Now need to summarise the MCMC output ...\n")
draws <- summarise_batch_anaphase_reversals(identifier,jobset_str_list,
                                     dt = dt, use_parallel=use_parallel)
cat("All done.\n")
