library(tidyverse)
library(sf)
#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate greater cirlce distance between:
#     (1) All clusters
#     (2) All samples
#---------------------------------------------------------------------------------------------------------------------------------------

#....................................................................
# Cluster Level
#....................................................................
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))
ge.nosf <- ge
sf::st_geometry(ge.nosf) <- NULL


gc <- ge.nosf %>%
  dplyr::select(c("longnum", "latnum")) %>%
  geosphere::distm(x =., fun = geosphere::distHaversine)

dim(gc)
length( unique(ge$dhsclust) )
colnames(gc) <- rownames(gc) <- ge$dhsclust


gc.long <- broom::tidy(as.dist( gc )) %>%
  dplyr::rename(gcdistance = distance)

saveRDS(gc.long, file = "data/distance_data/greater_circle_distance_forclusters.rds")



#....................................................................
# Sample Level
#....................................................................
drcsmpls <- sf::st_as_sf(readRDS("data/distance_data/drcsmpls_foruse.rds"))

smpl.gc <- geosphere::distm(x =sf::as_Spatial(drcsmpls), fun = geosphere::distHaversine)
colnames(smpl.gc) <- rownames(smpl.gc) <-  drcsmpls$id
smpl.gc.long <- broom::tidy(as.dist( smpl.gc )) %>%
  dplyr::rename(gcdistance = distance)

saveRDS(smpl.gc.long, file = "data/distance_data/greater_circle_distance_forsmpls.rds")

