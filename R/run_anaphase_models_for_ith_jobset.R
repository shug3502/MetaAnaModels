run_anaphase_models_for_ith_jobset <- function(i,jobset_str_list,K_list,fits_folder_str,identifier,
					       run_analysis,use_parallel,run_changept_anyway,dt,num_iter){
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
                                            fits_folder_str = fits_folder_str, dt=dt,
                                            plot_opt=0,num_iter=num_iter)
#t_ana_input_df <- changept_out %>%
#  spread_draws(t_ana[track]) %>% #spread draws is slow and uses lots of memory I think
#  summarise(t_ana=median(t_ana))

  t_ana_input_df <- as.data.frame(changept_out, pars = "t_ana") %>% 
    tidyr::gather(pair,t_ana) %>% 
    mutate(pair = as.numeric(stringr::str_extract(pair, "(\\d)+")))  %>% 
    group_by(pair) %>% 
    summarise(t_ana=median(t_ana))

out <- fit_anaphase_reversals_model(jobset_str_list[i],t_ana_input_df, K=K_list[i],
                                           identifier = identifier,
                                           run_analysis = run_analysis,
                                           use_parallel = use_parallel,
                                           fits_folder_str = fits_folder_str,
					   plot_opt=0, dt=dt,
					   num_iter=num_iter,
					   stan_file = here::here('src/stan_files/anaphase_reversals_hierarchical_bwd_sample.stan') 
)
cat("Done ...\n")
},
error = function (err) { print(err)
cat("Continuing but analysis not done for that job ... \n")
})
return(out)
}

