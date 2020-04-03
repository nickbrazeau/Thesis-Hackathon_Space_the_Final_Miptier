library(tidyverse)
library(sf)
#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this section is to calculate greater cirlce distance between clusters
#---------------------------------------------------------------------------------------------------------------------------------------
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

gc <- sf::st_distance(ge, which = "Great Circle")


dim(gc)
length( unique(ge$dhsclust) )
colnames(gc) <- rownames(gc) <- ge$dhsclust

gc.long <- broom::tidy(as.dist( gc )) %>%
  dplyr::rename(gcdistance = distance)

saveRDS(gc.long, file = "data/distance_data/greater_circle_distance_forclusters.rds")

