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

wtr.connected <- readRDS("data/derived_data/waterway_connected.rds")


ge2line <- geosphere::dist2Line(p = sf::as_Spatial(ge),
                                line = sf::as_Spatial(wtr.connected),
                                distfun = geosphere::distHaversine)

saveRDS(object = ge2line, file = "data/derived_data/ge2line.rds")
