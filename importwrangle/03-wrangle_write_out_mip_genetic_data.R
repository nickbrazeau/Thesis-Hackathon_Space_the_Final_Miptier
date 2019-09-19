#-----------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to write out the big barcode data
#-----------------------------------------------------------------------------------------------------------------------------------------
library(vcfR)
library(MIPanalyzer)
library(rplasmodium)

#.................................................
# Read in and force to monoclonal like Big Barcode MS did
# Make compatible with pf3d7
#.................................................

mipbivcfR <- MIPanalyzerbi2vcfR(input = readRDS(file = "data/raw_data/bigbarcode_genetic_data/biallelic_processed.rds"),
                                cutoff = 0.5)

liftover <- tibble(chrom = rplasmodium::chromnames(genome = "pf3d7")[1:14],
                   chr = paste0("chr", seq(1:14)))
CHROMbb <- left_join(tibble(chr = vcfR::getCHROM(mipbivcfR)), liftover)
mipbivcfR@fix[,1] <- unlist(CHROMbb[,2])


# write out
vcfR::write.vcf(x=mipbivcfR, file = "data/raw_data/bigbarcode_genetic_data/biallelic_processed.vcf.gz")

#....................
# Write out BB sample metadata
#....................
smplmtdt <- bigbarcode$samples
saveRDS(smplmtdt, file = "data/derived_data/sample_metadata.rds")
