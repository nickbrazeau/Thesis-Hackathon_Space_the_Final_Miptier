#-----------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to tidy the genetic data outputs that we
# created and subset the VCFs and IBD/S metrics
#-----------------------------------------------------------------------------------------------------------------------------------------
library(vcfR)
library(vcfRmanip)
library(tidyverse)

mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(country == "DRC")

# subset to samples with non-missing GPS data
dt <- readRDS("data/derived_data/DHS_qPCR_allkids_geo.rds")

mtdt <- mtdt %>%
  dplyr::filter(hv001 %in% dt$hv001)


#............................................................
# Subset vcf
#............................................................
vcf <- vcfR::read.vcfR(file = "data/derived_data/bigbarcode_genetic_data/polarized_biallelic_processed.wsaf.vcf.bgz")
vcf.DRC <- vcfRmanip::select_samples(vcf, smplvctr = unlist(mtdt$id))
saveRDS(object = vcf.DRC, file = "data/derived_data/bigbarcode_genetic_data/polarized_biallelic_processed_wsaf_DRC_VCFR.rds")

#............................................................
# Subset and Tidy IBD/S matrices
# from MIPANALYZER
#............................................................
DRCmp <- MIPanalyzer::vcf2mipanalyzer_biallelic(file = "data/derived_data/bigbarcode_genetic_data/polarized_biallelic_processed.wsaf.vcf.bgz")
# filter to DRC samples
drc.smpls <- DRCmp$samples$SAMPLE_ID %in% mtdt$id

DRCmp <- DRCmp %>%
  MIPanalyzer::filter_samples(., sample_filter = drc.smpls, description = "Subset to DRC Samples")

# get ibD
DRCmp.ibD <- MIPanalyzer::inbreeding_mle(DRCmp, f = seq(0.001, 0.999, 0.001),
                                         ignore_het = FALSE)
diag(DRCmp.ibD$mle) <- 1
colnames(DRCmp.ibD$mle) <- rownames(DRCmp.ibD$mle) <- mtdt$id

DRCmp.ibD.long <- broom::tidy(as.dist(t(DRCmp.ibD$mle))) %>%  # note, Bob returns upper triangle
  magrittr::set_colnames(c("smpl1", "smpl2", "malecotf"))

saveRDS(DRCmp.ibD.long, "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD_polarized_biallelic_processed.long.rds")



# get ibS
DRCmp.ibS <- MIPanalyzer::get_IBS_distance(x = DRCmp, ignore_het = F)

diag(DRCmp.ibS) <- 1
colnames(DRCmp.ibS) <- rownames(DRCmp.ibS) <- mtdt$id

DRCmp.ibS.long <- broom::tidy(as.dist(t(DRCmp.ibS))) %>%  # note, Bob returns upper triangle
  magrittr::set_colnames(c("smpl1", "smpl2", "hammings"))

saveRDS(DRCmp.ibS.long, "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS_polarized_biallelic_processed.long.rds")

#............................................................
# Save out the samples we will need for the sample level
# distance matrices
#............................................................
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::select(c("hv001", "geometry"))

drcsmpls <- mtdt %>%
  dplyr::select(c("id", "hv001")) %>%
  dplyr::left_join(., ge, by = "hv001")

saveRDS(drcsmpls, "data/distance_data/drcsmpls_foruse.rds")

