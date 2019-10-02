#................................................................................................................................
# Purpose of this script is to clean up distance measures and combine them
#................................................................................................................................
library(tidyverse)

#.............................
# Cluster Level
#.............................
gcclust <- readRDS("data/distance_data/greater_circle_distance_forclusters.rds") %>%
  dplyr::mutate(item1 = as.integer(item1),
                item2 = as.integer(item2))

roadsclust <- readRDS("data/distance_data/clstr_road_distmeters_long.rds")
# TODO water

distancesmatrix_cluster <- dplyr::inner_join(gcclust, roadsclust, by = c("item1", "item2")) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y", "gcdistance", "roaddistance"))




#................
# Out
#................
saveRDS(object = distancesmatrix_cluster, file = "data/distance_data/distancematrix_bycluster.rds")

