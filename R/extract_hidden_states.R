extract_hidden_states <- function(estimate,jobset_str,fits_folder_str,identifier) {
#give this function either a stanfit object, estimate, or the identifier to the job so can load in the fit

##########

if (is.null(estimate)) {
  cat("No estimate provided. Loading from saved fits objects...\n")
  job_id = stringr::str_split(jobset_str,"kittracking")[[1]][2]
  identifier = paste(job_id %>%
                       stringr::str_replace_all("\\.",""),identifier,sep="")
  path_to_fit <- file.path(fits_folder_str,paste('anaphase_reversals_hierarchical_',identifier,'.rds',sep=''))
  stopifnot(file.exists(path_to_fit))
  estimate = readRDS(path_to_fit)
  cat("Successfully loaded\n")
}

sigma_sim_matrix <- rstan::extract(estimate,pars="sigma_sim")$sigma_sim
#put in format that analysis code expects. That is a list of matrices, rather than multi dim array
niter <- dim(sigma_sim_matrix)[1]
sigma_sim <- purrr::map(1:niter,function(x) sigma_sim_matrix[x,,])
return(sigma_sim)
}
