#! /bin/bash

ROOT=/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/ # root directory for project (non-scratch)
SNAKE=/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/analyses/cluster_inbreeding
NODES=1028 # max number of cluster nodes
WAIT=30 # lag for system

snakemake \
	--snakefile $SNAKE/run_snake_spat_grad_descent.py \
	--configfile $SNAKE/config_clust.yaml \
	--directory $ROOT \
	--printshellcmds \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--cluster $SNAKE/launch.py \
	-j $NODES \
	--nolock \
#	--dryrun -p

# voroni prov tesselations
snakemake \
	--snakefile $SNAKE/run_snake_spat_grad_descent.py \
	--configfile $SNAKE/config_prov.yaml \
	--directory $ROOT \
	--printshellcmds \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--cluster $SNAKE/launch.py \
	-j $NODES \
	--nolock \
#	--dryrun -p
