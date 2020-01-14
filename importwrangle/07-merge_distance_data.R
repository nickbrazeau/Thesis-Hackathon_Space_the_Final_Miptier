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
riverclust <- readRDS("data/distance_data/river_distance_forclusters.rds") %>%
  dplyr::select(c("hv001.x", "hv001.y", "riverdist"))
riverclust$riverdist <- as.vector(unlist(riverclust$riverdist))

# several clusters appear to be snapped to the same river edge
# going to offset them by 5,000 meters which is near the min
riverclust <- riverclust %>%
  dplyr::mutate(riverdist = ifelse(riverdist == 0, 5e3, riverdist) )



distancesmatrix_cluster <- dplyr::inner_join(gcclust, roadsclust, by = c("item1", "item2")) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y", "gcdistance", "roaddistance"))

distancesmatrix_cluster <- long_distance_matrix_join(distancesmatrix_cluster, riverclust, by = c("hv001.x", "hv001.y"))
distancesmatrix_cluster$riverdist <- as.vector(unlist(distancesmatrix_cluster$riverdist))

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

roadrprov <- readRDS("data/distance_data/prov_road_distmeters_long.rds") %>%
  dplyr::mutate(item1 = as.integer(item1),
                item2 = as.integer(item2))


distancesmatrix_prov <- dplyr::inner_join(gcprov, roadrprov, by = c("item1", "item2")) %>%
  magrittr::set_colnames(c("item1", "item2", "gcdistance", "roaddistance"))



# liftover names
adm.dict <- tibble(item1 = 1:26,
                   adm1name.x = DRCprov$adm1name)
distancesmatrix_prov <- dplyr::left_join(distancesmatrix_prov, adm.dict, by = "item1")

adm.dict <- tibble(item2 = 1:26,
                   adm1name.y = DRCprov$adm1name)
distancesmatrix_prov <- dplyr::left_join(distancesmatrix_prov, adm.dict, by = "item2") %>%
  dplyr::select(-c("item1", "item2"))




riverprov <- readRDS("data/distance_data/river_distance_forprovinces.rds") %>%
  dplyr::select(c("from", "to", "riverdist")) %>%
  dplyr::rename(adm1name.x = from,
                adm1name.y = to)

distancesmatrix_prov <- long_distance_matrix_join(distancesmatrix_prov, riverprov,
                                                  by = c("adm1name.x", "adm1name.y"))

distancesmatrix_prov <- distancesmatrix_prov %>%
  dplyr::select(c("adm1name.x", "adm1name.y", "gcdistance", "roaddistance", "riverdist"))




#................
# Out
#................
saveRDS(object = distancesmatrix_prov, file = "data/distance_data/distancematrix_byprovince.rds")



