#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=49512
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier"); source("analysis/04-BayesSPlambda_backend.R")'
