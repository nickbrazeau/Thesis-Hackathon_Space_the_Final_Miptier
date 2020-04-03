#########################################################################
# Purpose: To calculate a distance matrix taking into account
#          flights. Going to "snap" a cluster to the nearest airport
#          and then weight their distance by the flow between airports
#
# Author: Nicholas F. Brazeau
#
# Date: April 01 2020
#########################################################################
library(tidyverse)
library(sf)
DRC <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds"))

#..................
# First pass, just assume you
# teleport --> will download aiprots from
# OSM Humanitarian
#..................
# https://data.humdata.org/dataset/ourairports-cod
airports <- readr::read_csv("data/raw_data/flight_data/hotosm_cd-airports.csv") %>%
  dplyr::filter(type %in% c("large_airport", "medium_airport")) %>%
  dplyr::select(c("name", "longitude_deg", "latitude_deg"))
airports <- sf::st_as_sf(airports, coords = c("longitude_deg", "latitude_deg"), crs = 4326)
ggplot() +
  geom_sf(data = DRC) +
  geom_sf(data = airports, color = "red")


#..................
# pull in clusters
#..................
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))
ge.osrm <- ge %>%
  dplyr::select("dhsclust")
rownames(ge.osrm) <- ge.osrm$dhsclust

#..................................
# SPIN up Docker and start server
# then add these options
#..................................
# https://github.com/rCarto/osrm
remotes::install_github("rCarto/osrm")
library(osrm)
options(osrm.server = "http://0.0.0.0:5000/", osrm.profile = "driving")

find_nearest_airport <- function(clst){
  #..............................................................
  # calculate minimum distance to an aiport
  #..............................................................
  ret <- min(osrm::osrmTable(src = clst, dst = airports,
                             measure = "distance")$distances) # in meters
  names(ret) <- clst$dhsclust
  return(ret)
}

ge.osrm <- split(ge.osrm, 1:nrow(ge.osrm))
ge.airportdist <- lapply(ge.osrm, find_nearest_airport)


ge.airportdist <- tibble::tibble(hv001 = unlist(lapply(ge.airportdist, names)),
                                 airportdist = unlist(ge.airportdist)) %>%
  dplyr::mutate(hv001 = as.numeric(hv001))


#..........
# Note, cluter 469 cannot be resolved with osrm
# do have one sample from that cluster...
# however cluster 313 is nearby, approximately 19763.25
# just going to add this approximate greater circle distance (2e4)
# to all 313 road distances to get 469 road dists
#..........
ge.airportdist$airportdist[ge.airportdist$hv001 == "469"] <- ge.airportdist$airportdist[ge.airportdist$hv001 == "313"] + 2e4

#..............................................................
# expand out and add in airport distances before being
# "teleported"
#..............................................................
clst_airport_distmat <- tibble::as_tibble(t(combn(ge$dhsclust, m = 2))) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y"))

colnames(ge.airportdist)[1] <- "hv001.x"
clst_airport_distmat <- dplyr::left_join(clst_airport_distmat, ge.airportdist, by = "hv001.x")
colnames(ge.airportdist)[1] <- "hv001.y"
clst_airport_distmat <- dplyr::left_join(clst_airport_distmat, ge.airportdist, by = "hv001.y")

#..................
# add and then drop
#..................
clst_airport_distmat <- clst_airport_distmat %>%
  dplyr::mutate(airport_distance = airportdist.x + airportdist.y) %>%
  dplyr::select(-c("airportdist.x", "airportdist.y"))

#..................
# save out
#..................
saveRDS(clst_airport_distmat, file = "data/distance_data/clstr_airport_distmeters_long.rds")







#..............................................................
# TATUM DATA
# read in flight data
#..............................................................
# arpts <- readr::read_csv("data/raw_data/flight_data/Airports_DRC.csv")
# flghts <- readr::read_csv("data/raw_data/flight_data/Air travel data in thousand_DRC.csv") %>%
#   magrittr::set_colnames(tolower(colnames(.))) %>%
#   dplyr::filter(origin_iso3 == "COD" & dest_iso3 == "COD")  # subset to DRC only
#
# # all flights together
# all( sort(unique(flghts$dest_airport)) %in% sort(unique(flghts$origin_airport)) )
# all( sort(unique(flghts$origin_airport)) %in% sort(unique(flghts$dest_airport)) )
# length(unique(flghts$origin_airport))
#
#
# length(unique(flghts$origin_airport))
#
#
# linker <- t(combn(unique(flghts$origin_airport), 2)) %>%
#   tibble::as_tibble(.) %>%
#   magrittr::set_colnames(c("arpt1", "arpt2")) %>%
#   dplyr::mutate(groupcode = 1:nrow(.))
# linker2 <- linker
# colnames(linker2) <- c("arpt2", "arpt1")
#
# t(combn(unique(flghts$origin_airport), 2))
