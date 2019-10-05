#............................................................................................................
# Purpose of this script is to run on cluster and get
# the nearest point along a line segment with respect to each dhsclust
#............................................................................................................
library(sf)
library(tidyverse)
library(geosphere)

ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

drc.rivers <-  sf::st_read("data/raw_data/drc_rivers_simplified/drc_rivers_simplified_postqgis.shp")

ge2line <- geosphere::dist2Line(p = sf::as_Spatial(ge),
                                line = sf::as_Spatial(drc.rivers),
                                distfun = geosphere::distHaversine)

saveRDS(object = ge2line, file = "data/derived_data/ge2line.rds")
