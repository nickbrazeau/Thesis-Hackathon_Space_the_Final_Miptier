#!/usr/bin/env bash
orig=~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data
path=~/Documents/GitHub/Space_the_Final_Miptier/analyses/admixture # user can define path

mkdir -p $path/gendat
gunzip -c $orig/bigbarcode_genetic_data/mipbivcfR.DRC.vcf.gz | bgzip > $path/gendat/mipbivcfR.DRC.vcf.bgz #written out as gzip compress, need bgzip
bcftools index $path/gendat/mipbivcfR.DRC.vcf.bgz
bcftools concat $path/gendat/mipbivcfR.DRC.vcf.bgz | bcftools sort | vcfdo wsaf | vcfdo dist > $path/mipbivcfR_DRC_dab.tab.txt
bcftools roh -G30 --estimate-AF - --rec-rate 7.4e-7 $path/gendat/mipbivcfR.DRC.vcf.bgz > $path/roh_estimates_DRC.bcftoolsout.txt
