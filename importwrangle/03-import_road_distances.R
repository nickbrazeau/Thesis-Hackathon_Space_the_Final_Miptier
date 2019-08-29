#................................................................................................................................
# Purpose of this script is to use OSRM to calculate
# road network distances between clusters
#................................................................................................................................
library(tidyverse)
library(sf)

# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

#..................................
# SPIN up Docker and start server
# then add these options
#..................................
# https://github.com/rCarto/osrm
remotes::install_github("rCarto/osrm")
library(osrm)
options(osrm.server = "http://0.0.0.0:5000/", osrm.profile = "driving")

ge.osrm <- ge %>%
  dplyr::select("dhsclust")
rownames(ge.osrm) <- ge.osrm$dhsclust


# now access API
clstrtrvldists <- osrm::osrmTable(loc = ge.osrm,
                                  measure = "distance") # in meters

#..........
# Note, cluter 469 cannot be resolved with osrm
# do have one sample from that cluster...
# however cluster 313 is nearby, approximately 19763.25
# just going to add this approximate greater circle distance (2e4)
# to all 313 road distances to get 469 road dists
#..........
clstrtrvldists.long <- broom::tidy(as.dist( clstrtrvldists$distances ))

clstr313 <- clstrtrvldists.long %>%
  dplyr::filter(item1 == 313 | item2 == 313) %>%
  dplyr::mutate(newitem1 = as.numeric( ifelse(item1 == 313, item1, item2) ),
                newitem2 = as.numeric( ifelse(item2 == 313, item1, item2))
                ) %>%
  dplyr::arrange(newitem2)

clstr313$distance <- clstr313$distance + 2e4

# don't rearrange here
clstr469 <- clstrtrvldists.long %>%
  dplyr::filter(item1 == 469 | item2 == 469) %>%
  dplyr::mutate(newitem2 = as.numeric( ifelse(item1 == 469, item2, item1) ) )


clstr469 <- left_join(clstr469, clstr313, by = "newitem2")

#..................
# now put dist back
#..................
clstrtrvldists.long$distance[is.na(clstrtrvldists.long$distance)] <- clstr469$distance.y


saveRDS(clstrtrvldists.long, file = "data/distance_data/clstr_travel_distmeters_long.rds")


