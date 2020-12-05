#-----------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate IBD and IBS
# from the Big Barcode Data
# Am also going to do some liftover to make a VCF (if I do downstream items)
# and make chrom names more compatible
#-----------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(vcfR)
remotes::install_github("nickbrazeau/vcfRmanip")
library(vcfRmanip)
remotes::install_github("mrc-ide/MIPanalyzer", ref = "master")
library(MIPanalyzer)
remotes::install_github("andrewparkermorgan/rplasmodium")
library(rplasmodium)

#.................................................
# Read in and force to monoclonal like Big Barcode MS did
# Make compatible with pf3d7
#.................................................
mipbi.mip <- readRDS(file = "data/raw_data/bigbarcode_genetic_data/biallelic_processed.rds")
liftover <- tibble(chrompf = rplasmodium::chromnames(genome = "pf3d7")[1:14],
                   CHROM = paste0("chr", seq(1:14)))
mipbi.mip$loci <- dplyr::left_join(mipbi.mip$loci, liftover, by = "CHROM")
mipbi.mip$loci <- mipbi.mip$loci %>%
  dplyr::select(-c("CHROM")) %>%
  dplyr::rename(CHROM = chrompf)

#..............................................................
# Subset to DRC Samples
#..............................................................
mtdt <- readRDS("data/derived_data/sample_metadata.rds")

# filter to DRC samples
# need boolean vector
drc.smpls <- mipbi.mip$samples$ID %in% mtdt$name

mipbi.mip.DRC <- mipbi.mip %>%
  MIPanalyzer::filter_samples(., sample_filter = drc.smpls, description = "Subset to DRC Samples")


#..............................................................
# Calculate IBD
#..............................................................
DRCmp.ibD <- MIPanalyzer::inbreeding_mle(mipbi.mip.DRC,
                                         f = seq(0.001, 0.999, 0.001),
                                         ignore_het = FALSE)
diag(DRCmp.ibD$mle) <- 1
colnames(DRCmp.ibD$mle) <- rownames(DRCmp.ibD$mle) <- mtdt$name

DRCmp.ibD.long <- broom::tidy(as.dist(t(DRCmp.ibD$mle))) %>%  # note, Bob returns upper triangle
  magrittr::set_colnames(c("smpl1", "smpl2", "malecotf"))

dir.create("data/derived_data/bigbarcode_genetic_data/", recursive = T)
saveRDS(DRCmp.ibD.long, "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.rds")

#..............................................................
# Calculate IBS
#..............................................................
DRCmp.ibS <- MIPanalyzer::get_IBS_distance(x = mipbi.mip.DRC,
                                           ignore_het = F)

diag(DRCmp.ibS) <- 1
colnames(DRCmp.ibS) <- rownames(DRCmp.ibS) <- mtdt$name

DRCmp.ibS.long <- broom::tidy(as.dist(t(DRCmp.ibS))) %>%  # note, Bob returns upper triangle
  magrittr::set_colnames(c("smpl1", "smpl2", "hammings"))

saveRDS(DRCmp.ibS.long, "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS.long.rds")

