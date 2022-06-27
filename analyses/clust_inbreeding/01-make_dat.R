#########################################################################
# Purpose: Wrangle data for cluster inbreeding grad descent function
#
# Date: April 02 2020
#########################################################################
source("R/pairwise_helpers.R")
library(tidyverse)

#............................................................
#### data for all  #####
#...........................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")

#......................
# distance data
#......................
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

airdist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "airportdistance")) %>%
  dplyr::mutate(airportdistance = ifelse(hv001.x == hv001.y, 0, airportdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))

#......................
# save out
#......................
dir.create("data/derived_data/allsmpls_clst_inbreeding_dat")
saveRDS(as.data.frame(gcdist_gens),
        file = "data/derived_data/allsmpls_clst_inbreeding_dat/gcdist_gens.RDS")
saveRDS(as.data.frame(roaddist_gens),
        file = "data/derived_data/allsmpls_clst_inbreeding_dat/roaddist_gens.RDS")
saveRDS(as.data.frame(airdist_gens),
        file = "data/derived_data/allsmpls_clst_inbreeding_dat/airdist_gens.RDS")


#............................................................
#### data for monoclonals  #####
#...........................................................
#......................
# subset to monoclonals
#......................
monoclonals <- readRDS("data/derived_data/RMCL_results/monoclonals.rds")

ibD.monoclonals <- ibD %>%
  dplyr::filter(barcode.x %in% monoclonals$barcode & barcode.y %in% monoclonals$barcode)

#......................
# distance data
#......................
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster) %>%
  dplyr::mutate(hv001.x = as.numeric(as.character(hv001.x)),
                hv001.y = as.numeric(as.character(hv001.y))) # factor liftover

gcdist_gens_monoclonals <- ibD.monoclonals %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "gcdistance")) %>%
  dplyr::mutate(gcdistance = ifelse(hv001.x == hv001.y, 0, gcdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


roaddist_gens_monoclonals <- ibD.monoclonals %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "roaddistance")) %>%
  dplyr::mutate(roaddistance = ifelse(hv001.x == hv001.y, 0, roaddistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))

airdist_gens_monoclonals <- ibD.monoclonals %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "airportdistance")) %>%
  dplyr::mutate(airportdistance = ifelse(hv001.x == hv001.y, 0, airportdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))

#......................
# save out
#......................
dir.create("data/derived_data/coione_clst_inbreeding_dat")
saveRDS(as.data.frame(gcdist_gens_monoclonals),
        file = "data/derived_data/coione_clst_inbreeding_dat/gcdist_gens.RDS")
saveRDS(as.data.frame(roaddist_gens_monoclonals),
        file = "data/derived_data/coione_clst_inbreeding_dat/roaddist_gens.RDS")
saveRDS(as.data.frame(airdist_gens_monoclonals),
        file = "data/derived_data/coione_clst_inbreeding_dat/airdist_gens.RDS")





