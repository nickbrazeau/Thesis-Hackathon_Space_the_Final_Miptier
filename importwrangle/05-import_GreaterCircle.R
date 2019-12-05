library(tidyverse)
library(sf)
#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this section is to calculate greater cirlce distance between clusters
#---------------------------------------------------------------------------------------------------------------------------------------

ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))
ge.nosf <- ge
sf::st_geometry(ge.nosf) <- NULL


gc <- ge.nosf %>%
  dplyr::select(c("longnum", "latnum")) %>%
  geosphere::distm(x =., fun = geosphere::distGeo)

dim(gc)
length( unique(ge$dhsclust) )
colnames(gc) <- rownames(gc) <- ge$dhsclust


gc.long <- broom::tidy(as.dist( gc )) %>%
  dplyr::rename(gcdistance = distance)

saveRDS(gc.long, file = "data/distance_data/greater_circle_distance_forclusters.rds")

#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this section is to calculate greater cirlce distance between province
#---------------------------------------------------------------------------------------------------------------------------------------
drcpov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")) %>%
  dplyr::mutate(provcentroid = sf::st_centroid(geometry))
drcpov.nosf <- drcpov
sf::st_geometry(drcpov.nosf) <- NULL

gc <- drcpov.nosf %>%
  dplyr::mutate(longnum = sf::st_coordinates(provcentroid)[,1],
                latnum = sf::st_coordinates(provcentroid)[,2]) %>%
  dplyr::select(c("longnum", "latnum")) %>%
  geosphere::distm(x =., fun = geosphere::distGeo)

dim(gc)
length( unique(drcpov.nosf$adm1name) )
colnames(gc) <- rownames(gc) <- drcpov.nosf$adm1name


gc.long <- broom::tidy(as.dist( gc )) %>%
  dplyr::rename(gcdistance = distance)




saveRDS(gc.long, file = "data/distance_data/greater_circle_distance_forprovinces.rds")


