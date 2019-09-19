
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



#............................................................
# Subset vcf
#............................................................
vcf <- vcfR::read.vcfR(file = "data/derived_data/bigbarcode_genetic_data/polarized_biallelic_processed.wsaf.vcf.bgz")
vcf.DRC <- vcfRmanip::select_samples(vcf, smplvctr = unlist(mtdt$id))
saveRDS(object = vcf.DRC, file = "data/derived_data/bigbarcode_genetic_data/polarized_biallelic_processed_wsaf_DRC_VCFR.rds")

#............................................................
# Subset and Tidy IBD/S matrices
#............................................................

ibd <- readr::read_tsv("data/derived_data/bigbarcode_genetic_data/IBD_polarized_biallelic_processed.long.tab.txt", col_names = F) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "nsites", "malecotf")) %>%
  dplyr::filter(smpl1 %in% mtdt$id & smpl2 %in% mtdt$id)

saveRDS(ibd, "data/derived_data/bigbarcode_genetic_data/DRCIBD_polarized_biallelic_processed.long.rds")

ibs <- read.table("data/derived_data/bigbarcode_genetic_data/IBS_polarized_biallelic_processed.long.tab.txt",
                  header = T, stringsAsFactors = F)
rownames(ibs) <- unlist(ibs[,1])
ibs <- ibs[,-1]
ibs <- broom::tidy( as.dist(ibs) ) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "hammingibs"))  %>%
  dplyr::filter(smpl1 %in% mtdt$id & smpl2 %in% mtdt$id)

saveRDS(ibd, "data/derived_data/bigbarcode_genetic_data/DRCIBS_polarized_biallelic_processed.long.rds")




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

