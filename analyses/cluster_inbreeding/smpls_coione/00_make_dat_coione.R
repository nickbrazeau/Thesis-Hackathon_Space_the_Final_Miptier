#########################################################################
# Purpose: Wrangle data for cluster inbreeding grad descent function
#
# Date: April 02 2020
#########################################################################
source("R/pairwise_helpers.R")
library(tidyverse)

#......................
# wrangle and take to long, diagonals to 0 for distances
#......................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")

#......................
# subset to monoclonals
#......................
coi <- readRDS("data/raw_data/RMCL_results/non_summariased_cois.rds") %>%
  dplyr::filter(region_denom == 1) %>% # just subset to DRc
  dplyr::filter(region == "DRC") %>%
  dplyr::filter(gt == 0.1) %>% # same as in Big Barcode Manuscript
  dplyr::mutate(barcode = ifelse(nchar(name) == 9, paste(strsplit(name, split = "")[[1]][5:nchar(name)], collapse = ""),
                                 ifelse(nchar(name) == 8, paste(strsplit(name, split = "")[[1]][4:nchar(name)], collapse = ""),
                                        ifelse(nchar(name) == 6, paste(strsplit(name, split = "")[[1]][2:nchar(name)], collapse = ""),
                                               name))))

# metadata
drcsmpls <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("barcode", "hv001", "longnum", "latnum"))

# bring together
coi_ge <- dplyr::left_join(coi, drcsmpls, by = "barcode")
monoclonals <- coi_ge %>% dplyr::filter(median == 1)
dir.create("data/derived_data/RMCL_results/", recursive = T)
saveRDS(monoclonals, "data/derived_data/RMCL_results/monoclonals.rds")

ibD.monoclonals <- ibD %>%
  dplyr::filter(barcode.x %in% monoclonals$barcode & barcode.y %in% monoclonals$barcode)

#......................
# distance data
#......................
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster) %>%
  dplyr::mutate(hv001.x = as.numeric(as.character(hv001.x)),
                hv001.y = as.numeric(as.character(hv001.y))) # factor liftover

gcdist_gens <- ibD.monoclonals %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "gcdistance")) %>%
  dplyr::mutate(gcdistance = ifelse(hv001.x == hv001.y, 0, gcdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


roaddist_gens <- ibD.monoclonals %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "roaddistance")) %>%
  dplyr::mutate(roaddistance = ifelse(hv001.x == hv001.y, 0, roaddistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))

airdist_gens <- ibD.monoclonals %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "airportdistance")) %>%
  dplyr::mutate(airportdistance = ifelse(hv001.x == hv001.y, 0, airportdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))




#......................
# save out
#......................
dir.create("data/derived_data/coione_clst_inbreeding_dat")
saveRDS(as.data.frame(gcdist_gens),
        file = "data/derived_data/coione_clst_inbreeding_dat/gcdist_gens.RDS")
saveRDS(as.data.frame(roaddist_gens),
        file = "data/derived_data/coione_clst_inbreeding_dat/roaddist_gens.RDS")
saveRDS(as.data.frame(airdist_gens),
        file = "data/derived_data/coione_clst_inbreeding_dat/airdist_gens.RDS")

