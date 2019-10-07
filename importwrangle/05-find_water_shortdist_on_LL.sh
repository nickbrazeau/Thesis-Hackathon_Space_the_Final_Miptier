#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=49512
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

R -e "setwd('/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier'); source('importwrangle/05-05-slurmfind_short_river_dist.R')"
