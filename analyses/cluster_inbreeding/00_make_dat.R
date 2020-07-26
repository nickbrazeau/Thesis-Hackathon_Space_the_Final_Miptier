#########################################################################
# Purpose: Wrangle data for cluster inbreeding grad descent function
#
# Date: April 02 2020
#########################################################################
source("R/pairwise_helpers.R")
library(tidyverse)

#......................
# wrangle and take to long, diagonals to 0
#......................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster)

gcdist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "gcdistance")) %>%
  dplyr::mutate(gcdistance = ifelse(hv001.x == hv001.y, 0, gcdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


roaddist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "roaddistance")) %>%
  dplyr::mutate(roaddistance = ifelse(hv001.x == hv001.y, 0, roaddistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


riverdist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "riverdistance")) %>%
  dplyr::mutate(riverdistance = ifelse(hv001.x == hv001.y, 0, riverdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))



airplanedist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "airport_distance")) %>%
  dplyr::mutate(airport_distance = ifelse(hv001.x == hv001.y, 0, airport_distance),) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))

#......................
# save out
#......................
dir.create("data/derived_data/clst_inbreeding_dat")
saveRDS(as.data.frame(gcdist_gens),
        file = "data/derived_data/clst_inbreeding_dat/gcdist_gens.RDS")
saveRDS(as.data.frame(roaddist_gens),
        file = "data/derived_data/clst_inbreeding_dat/roaddist_gens.RDS")
saveRDS(as.data.frame(riverdist_gens),
        file = "data/derived_data/clst_inbreeding_dat/riverdist_gens.RDS")
saveRDS(as.data.frame(airplanedist_gens),
        file = "data/derived_data/clst_inbreeding_dat/airplanedist_gens.RDS")



