#!/bin/bash 

sbatch -o output/out0818_missing.log --export=ALL,path_to_data='data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame' src/bash/run_mcmc_revision.sh
sbatch -o output/out0820_missing.log --export=ALL,path_to_data='data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame' src/bash/run_mcmc_revision.sh
sbatch -o output/out0929_missing.log --export=ALL,path_to_data='data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame' src/bash/run_mcmc_revision.sh
sbatch -o output/out0930_missing.log --export=ALL,path_to_data='data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame' src/bash/run_mcmc_revision.sh
