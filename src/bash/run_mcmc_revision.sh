#!/bin/bash
#SBATCH --job-name=MetaAnaModels_2s_revision                          #see pauls email or eg here https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts
#SBATCH --time=180:00:00
#SBATCH --mem-per-cpu=8g
#SBATCH --export=ALL
#SBATCH --ntasks=1                                        # Run on a single CPU

pwd; hostname; date
module load R
identifier="${USER}_v332"
echo "Rscript poor_oscillators_anaphase_analysis_part3.R ${path_to_data} ${identifier}"
Rscript poor_oscillators_anaphase_analysis_part3.R "${path_to_data}" "${identifier}" 
