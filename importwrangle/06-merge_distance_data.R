#................................................................................................................................
# Purpose of this script is to clean up distance measures and combine them
#................................................................................................................................
library(tidyverse)

#.............................
# Cluster Level
#.............................
gclust <- readRDS("data/distance_data/greater_circle_distance_forclusters.rds") %>%
  dplyr::mutate(item1 = as.integer(item1),
                item2 = as.integer(item2))

roadsclust <- readRDS("data/distance_data/clstr_road_distmeters_long.rds")
# TODO water

distancesmatrix_cluster <- dplyr::inner_join(gclust, roadsclust, by = c("item1", "item2"))


#.............................
# Sample Level
#.............................
gcsmpls <- readRDS("data/distance_data/greater_circle_distance_forsmpls.rds")
roadssmpls <- readRDS("data/distance_data/smpl_road_distmeters_long.rds")

# TODO water


distancesmatrix_smpl <- dplyr::inner_join(gcsmpls, roadssmpls, by = c("item1", "item2"))

#................
# Out
#................

saveRDS(object = distancesmatrix_cluster, file = "data/distance_data/distancematrix_bycluster.rds")
saveRDS(object = distancesmatrix_smpl, file = "data/distance_data/distancematrix_bysmpl.rds")

