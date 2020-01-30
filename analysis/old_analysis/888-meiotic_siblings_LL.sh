#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=128g
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier"); source("analysis/05-IBD_meioticsibs_backend.R")'
