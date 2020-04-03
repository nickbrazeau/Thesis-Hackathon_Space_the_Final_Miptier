#! /bin/bash

ROOT=/Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/ # root directory for project (non-scratch)
SNAKE=/Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/analysis/cluster_inbreeding
WAIT=30 # lag for system

snakemake \
	--snakefile $SNAKE/run_snake_params.py \
	--configfile $SNAKE/config_temp.yaml \
	--directory $ROOT \
	--printshellcmds \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--dryrun -p
