library(tidyverse)
library(sf)
library(raster)
library(shp2graph)
#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate greater cirlce distance between:
#     (1) All clusters
#     (2) All samples
#---------------------------------------------------------------------------------------------------------------------------------------
#...............................
# Waterways
#...............................
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
wtrlns <- sf::read_sf("data/raw_data/hotosm_waterway/hotosm_cod_waterways_grass_qgis_cleaned.shp")
# pull out network of longest connection
wtr.connected <- shp2graph::nt.connect(sf::as_Spatial(wtrlns))
proj4string(wtr.connected) <- "+proj=longlat +datum=WGS84 +no_defs"

ggplot() +
  geom_sf(data=DRCprov) +
  geom_sf(data = wtrlns, color = "blue") +
  geom_sf(ge, color = "red")

# looks reasonable

#....................................................................
# Cluster Level
#....................................................................
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))


#............................................
# Snap clusters to rivers and find distance
#............................................
points <- ge[,c("longnum", "latnum")]
sf::st_geometry(points) <- NULL

network.river.ge <- shp2graph::points2network(ntdata = wtr.connected,
                                              pointsxy = points,
                                              approach = 4,
                                              Detailed = T,
                                              longlat = T,
                                              ELComputed = T,
                                              ea.prop = rep(0,13)) #https://rdrr.io/cran/shp2graph/man/points2network.html

shp2graph::ptsinnt.view(ntdata= wtr.connected,
                        nodelist=network.river.ge[[1]],
                        pointsxy=points,
                        VElist=network.river.ge[[7]],
                        CoorespondIDs=network.river.ge[[3]])


igraph.river.ge <- shp2graph::nel2igraph(nodelist = network.river.ge[[1]],
                                         edgelist = network.river.ge[[2]],
                                         weight =   network.river.ge[[8]])

#................
# Get Distances
#................
dhspoints <- as.numeric(network.river.ge[[3]])
rivermat <- igraph::shortest.paths(igraph.river.ge, v = dhspoints, to = dhspoints)
rivermat <- rivermat[, match(dhspoints, dhspoints)]
colnames(rivermat) <- rownames(rivermat) <- ge$dhsclust

#................
# Out
#................
river_long <- broom::tidy(as.dist( rivermat )) %>%
  dplyr::rename(waterdistance = distance)

saveRDS(object = river_long, file = "data/distance_data/cluster_water_distance.rds")

