#!/bin/bash 

sbatch -o output/out0818.log src/bash/run_mcmc_MetaAnaModels.sh --export=path_to_data="data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame"
sbatch -o output/out0820.log src/bash/run_mcmc_MetaAnaModels.sh --export=path_to_data="data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame"
sbatch -o output/out0929.log src/bash/run_mcmc_MetaAnaModels.sh --export=path_to_data="data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame"
sbatch -o output/out0930.log src/bash/run_mcmc_MetaAnaModels.sh --export=path_to_data="data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame"
