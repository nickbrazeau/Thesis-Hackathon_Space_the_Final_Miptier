#!/usr/bin/env bash

vcfindir=/Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/raw_data/bigbarcode_genetic_data/
vcfoutdir=/Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/
FASTA=/Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/ancestral.fa
# apm downloaded from p_reich MS and imputed alleles
# /proj/ideel/julianog/users/apm/p_reichenowi
# cited in the bip mip paper

# housekeeping
gunzip -c $vcfindir/biallelic_processed.vcf.gz | bgzip >  $vcfoutdir/biallelic_processed.vcf.bgz #written out as gzip compress, need bgzip
bcftools index $vcfoutdir/biallelic_processed.vcf.bgz

# polarize
bcftools concat $vcfoutdir/biallelic_processed.vcf.bgz  | \
bcftools sort | \
vcfdo polarize -f $FASTA | bgzip >  $vcfoutdir/polarized_biallelic_processed.vcf.bgz

# wsaf and fws
bcftools concat $vcfoutdir/polarized_biallelic_processed.vcf.bgz  | \
vcfdo wsaf | bgzip >  $vcfoutdir/polarized_biallelic_processed.wsaf.vcf.bgz
