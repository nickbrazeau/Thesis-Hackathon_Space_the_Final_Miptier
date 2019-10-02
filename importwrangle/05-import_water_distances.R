#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate the water distanc between clusters
# Following this tutorial https://www.r-spatial.org/r/2019/09/26/spatial-networks.html
#---------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(tidygraph)
library(igraph)
library(shp2graph)
library(raster)

#............................................................................................................
# IMPORT DATA
#............................................................................................................
#...............................
# Ge import
#...............................
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")

#...............................
# Waterways
#...............................
wtrlns <- sf::read_sf("data/raw_data/hotosm_waterway/hotosm_cod_waterways_grass_qgis_cleaned.shp")
# pull out network of longest connection
wtr.connected <- shp2graph::nt.connect(sf::as_Spatial(wtrlns))
proj4string(wtr.connected) <- "+proj=longlat +datum=WGS84 +no_defs"
saveRDS(object = wtr.connected, file = "data/derived_data/waterway_connected.rds") # for cluster

# quick view
ggplot() +
  geom_sf(data=DRCprov, color = "#f0f0f0") +
  geom_sf(data = wtr.connected, color = "#3288bd") +
  geom_sf(dat = ge, color = "#f22121") +
  theme_minimal()



#............................................................................................................
# From Connect DRC Cluster Points to Riverways
#............................................................................................................
# find nearest line per point
# note, ran this on the cluster and return here
ge2line <- readRDS("data/derived_data/ge2line.rds")

ge.start <- data.frame( sf::st_coordinates(ge) )
ge.end <- data.frame(X = ge2line$lon, Y = ge2line$lat)


clusterlines <- data.frame(osm_id = NA, name = "dhsclustnode",
                           waterway = NA, geometry = rep(NA, nrow(ge.end)) )

for(i in 1:nrow(ge.start)){
  start <- sf::st_point(x = c(ge.start$X[i], ge.start$Y[i]))
  end <- sf::st_point(x = c(ge.end$X[i], ge.end$Y[i]))
  startend <- sf::st_combine(c(start, end))
  clusterlines$geometry[i] <- sf::st_cast(startend, "LINESTRING")
}

clusterlines <- sf::st_as_sf(clusterlines)


# now combine
wtr.connected <- sf::st_as_sf(wtr.connected) %>%
  dplyr::select(c("osm_id", "name", "waterway", "geometry"))

wtr.connected <- rbind.data.frame(wtr.connected, clusterlines)

#............................................................................................................
# Manipulate Shapes to Prepare for Network
#............................................................................................................


#...............................
# Edges
# Give each line a unique ID
#...............................
edges <- wtr.connected %>%
  mutate(edgeID = c(1:n()))


#...............................
# Nodes
# give each end and start of line a unique ID
#...............................
nodes <- edges %>%
  sf::st_coordinates() %>%
  tibble::as_tibble() %>%
  dplyr::rename(edgeID = L1) %>%
  dplyr::group_by(edgeID) %>%
  dplyr::slice(c(1, n())) %>% # removing duplicate edge IDs that will be future node IDs
  dplyr::ungroup() %>%
  dplyr::mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
  dplyr::mutate(xy = paste(.$X, .$Y)) %>%
  dplyr::mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
  dplyr::select(-xy)


#...............................
# Collect start and end nodes
# and index to the edges
#...............................
source_nodes <- nodes %>%
  dplyr::filter(start_end == 'start') %>%
  dplyr::pull(nodeID)

target_nodes <- nodes %>%
  dplyr::filter(start_end == 'end') %>%
  dplyr::pull(nodeID)

edges = edges %>%
  dplyr::mutate(from = source_nodes, to = target_nodes)


#...............................
# Collect final nodes
#...............................
nodes <- nodes %>%
  dplyr::distinct(nodeID, .keep_all = TRUE) %>%
  dplyr::select(-c(edgeID, start_end)) %>%
  sf::st_as_sf(coords = c('X', 'Y')) %>%
  sf::st_set_crs(st_crs(edges))



#.............................
# Label cluster nodes
#.............................
# use sf::st intersect or somehing to figure out write nodes for clusters

#............................................................................................................
# Make Network
#............................................................................................................





