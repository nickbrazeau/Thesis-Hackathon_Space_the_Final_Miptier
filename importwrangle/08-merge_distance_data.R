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

# roads
roadsclust <- readRDS("data/distance_data/clstr_road_distmeters_long.rds")

# planes
planeclust <- readRDS(file = "data/distance_data/clstr_airport_distmeters_long.rds")

# river offset
riverclust <- readRDS("data/distance_data/river_distance_forclusters.rds") %>%
  dplyr::select(c("hv001.x", "hv001.y", "riverdistance"))

# several clusters appear to be snapped to the same river edge
# going to offset them by 5,000 meters which is near the min
riverclust <- riverclust %>%
  dplyr::mutate(riverdistance = ifelse(riverdistance == 0, 5e3, riverdistance) )



distancesmatrix_cluster <- dplyr::inner_join(gcclust, roadsclust, by = c("item1", "item2")) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y", "gcdistance", "roaddistance"))

distancesmatrix_cluster <- long_distance_matrix_join(distancesmatrix_cluster, riverclust,
                                                     by = c("hv001.x", "hv001.y"))
distancesmatrix_cluster$riverdistance <- as.vector(unlist(distancesmatrix_cluster$riverdistance))

distancesmatrix_cluster <- long_distance_matrix_join(distancesmatrix_cluster, planeclust,
                                                     by = c("hv001.x", "hv001.y"))

#................
# Out
#................
saveRDS(object = distancesmatrix_cluster, file = "data/distance_data/distancematrix_bycluster.rds")
