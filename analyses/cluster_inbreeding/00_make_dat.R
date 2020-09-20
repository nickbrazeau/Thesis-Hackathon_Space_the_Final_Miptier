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
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster) %>%
  dplyr::mutate(hv001.x = as.numeric(as.character(hv001.x)),
                hv001.y = as.numeric(as.character(hv001.y))) # factor liftover

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

#......................
# wrangle migration provinces demes
#......................
node_pairs <- readRDS(file = "data/distance_data/vr_nodepairs_migrate_disance.rds")

migrate_gens <- ibD %>%
  dplyr::left_join(., node_pairs, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "PrdMIG")) %>%
  dplyr::mutate(roaddistance = ifelse(hv001.x == hv001.y, 0, PrdMIG)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))



#......................
# save out
#......................
dir.create("data/derived_data/clst_inbreeding_dat")
saveRDS(as.data.frame(gcdist_gens),
        file = "data/derived_data/clst_inbreeding_dat/gcdist_gens.RDS")
saveRDS(as.data.frame(roaddist_gens),
        file = "data/derived_data/clst_inbreeding_dat/roaddist_gens.RDS")
saveRDS(as.data.frame(migrate_gens),
        file = "data/derived_data/clst_inbreeding_dat/migrate_gens.RDS")


