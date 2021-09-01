#!/bin/bash
#SBATCH --job-name=MetaAnaModels_2s_hierarchical			  #see pauls email or eg here https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts
#SBATCH --time=180:00:00
#SBATCH --mem-per-cpu=32g
#SBATCH --export=ALL
#SBATCH --ntasks=1                   			  # Run on a single CPU

pwd; hostname; date
module load R
identifier="${USER}_MetaAnaModels_2s_v421"
echo "Rscript R/main ${identifier} ${path_to_data}"
Rscript R/main.R "${identifier}" "${path_to_data}"
