#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=256g
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier"); source("analysis/04-MCMC_BayesBetween_backend.R")'
Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier"); source("analysis/04-MCMC_BayesWithin_backend.R")'
