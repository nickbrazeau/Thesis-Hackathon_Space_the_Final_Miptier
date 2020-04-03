#########################################################################
# Purpose: Make data that will go into the cluster inbreeding map
#
# Author: Nicholas F. Brazeau
#
# Date: April 02 2020
#########################################################################
source("R/pairwise_helpers.R")
library(tidyverse)
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster)

gcdist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "gcdistance")) %>%
  dplyr::mutate(gcdistance = ifelse(hv001.x == hv001.y, 0, gcdistance),
                gcdistance = gcdistance / sd(gcdistance)
                ) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


roaddist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "roaddistance")) %>%
  dplyr::mutate(roaddistance = ifelse(hv001.x == hv001.y, 0, roaddistance),
                roaddistance = roaddistance / sd(roaddistance)
                ) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


#..............................................................
# save out
#..............................................................
dir.create("data/derived_data/clst_inbreeding_dat")
saveRDS(as.data.frame(gcdist_gens),
        file = "data/derived_data/clst_inbreeding_dat/gcdist_gens.RDS")
saveRDS(as.data.frame(roaddist_gens),
        file = "data/derived_data/clst_inbreeding_dat/roaddist_gens.RDS")



