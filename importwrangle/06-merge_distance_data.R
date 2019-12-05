#................................................................................................................................
# Purpose of this section is to clean up distance measures and combine them
# for CLUSTERS
#................................................................................................................................
library(tidyverse)
source("R/pairwise_helpers.R")
#.............................
# Cluster Level
#.............................
gcclust <- readRDS("data/distance_data/greater_circle_distance_forclusters.rds") %>%
  dplyr::mutate(item1 = as.integer(item1),
                item2 = as.integer(item2))

roadsclust <- readRDS("data/distance_data/clstr_road_distmeters_long.rds")

# river offset
km1 <- 1e3
units(km1) <- units::as_units("m")
riverclust <- readRDS("data/distance_data/river_distance_forclusters.rds") %>%
  dplyr::select(c("dhsclustfrom", "dhsclustto", "riverdist")) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y", "riverdist")) %>%
  dplyr::mutate(riverdist = riverdist + km1) # river distance offset factor of 1km

distancesmatrix_cluster <- dplyr::inner_join(gcclust, roadsclust, by = c("item1", "item2")) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y", "gcdistance", "roaddistance"))

distancesmatrix_cluster <- long_distance_matrix_join(distancesmatrix_cluster, riverclust, by = c("hv001.x", "hv001.y"))

#................
# Out
#................
saveRDS(object = distancesmatrix_cluster, file = "data/distance_data/distancematrix_bycluster.rds")


#................................................................................................................................
# Purpose of this section is to clean up distance measures and combine them
# for PROVINCE
#................................................................................................................................
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")
gcprov <- readRDS("data/distance_data/greater_circle_distance_forprovinces.rds") %>%
  dplyr::mutate(item1 = as.integer(item1),
                item2 = as.integer(item2))

roadrprov <- readRDS("data/distance_data/prov_road_distmeters_long.rds")

riverprov <- readRDS("data/distance_data/river_distance_forprovinces.rds") %>%
  dplyr::select(c("item1", "item2", "riverdist"))  %>%
  dplyr::mutate(riverdist = riverdist + km1) # river distance offset factor of 1km

distancesmatrix_prov <- dplyr::inner_join(gcprov, roadrprov, by = c("item1", "item2")) %>%
  magrittr::set_colnames(c("item1", "item2", "gcdistance", "roaddistance"))

distancesmatrix_prov <- long_distance_matrix_join(distancesmatrix_prov, riverprov,
                                                  by = c("item1", "item2"))

# liftover names
adm.dict <- tibble(item1 = 1:26,
                   adm1name.x = DRCprov$adm1name)
distancesmatrix_prov <- dplyr::left_join(distancesmatrix_prov, adm.dict, by = "item1")

adm.dict <- tibble(item2 = 1:26,
                   adm1name.y = DRCprov$adm1name)
distancesmatrix_prov <- dplyr::left_join(distancesmatrix_prov, adm.dict, by = "item2") %>%
  dplyr::select(-c("item1", "item2"))


#................
# Out
#................
saveRDS(object = distancesmatrix_prov, file = "data/distance_data/distancematrix_byprovince.rds")



