#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate the water distanc between clusters
# Following this tutorial https://www.r-spatial.org/r/2019/09/26/spatial-networks.html
# to get shorest distances along the river network
# Above allows us to make a tidy spatial network
#---------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(tidygraph)
library(igraph)
library(shp2graph)
library(raster)
source("R/themes.R")
source("R/pairwise_helpers.R")

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
wtrlns <- sf::read_sf("data/raw_data/hotosm_cod_waterways_lines_shp/hotosm_cod_waterways_lines.shp")
wtrlns <- wtrlns %>%
  dplyr::filter(waterway == "river" & !is.na(name)) %>%
  dplyr::select(c("name", "geometry"))

wtrpolys <- sf::read_sf("data/raw_data/hotosm_cod_waterways_polygons_shp/hotosm_cod_waterways_polygons.shp")
wtrpolys <- wtrpolys %>%
  dplyr::filter(name == "Lake Tanganyika") %>%
  dplyr::select(c("name", "geometry"))

#............................................................................................................
# write out and manipulate in qgis
#............................................................................................................

drc.rivers <- rbind.data.frame(wtrlns, wtrpolys)

dir.create("data/raw_data/drc_rivers_simplified/")
sf::write_sf(obj = drc.rivers,
             dsn = "data/raw_data/drc_rivers_simplified/drc_rivers_init.shp") # for cluster

# read back in
drc.rivers <- sf::st_read("data/raw_data/drc_rivers_simplified/drc_rivers_simplified_postqgis.shp")


#............................................................................................................
# Manipulate Shapes to Prepare for Network
#............................................................................................................
#...............................
# Edges
# Give each line a unique ID
#...............................
edges <- drc.rivers %>%
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



#............................................................................................................
# Make (and plot) Network
#............................................................................................................
rivernetwork <- tidygraph::tbl_graph(nodes = nodes,
                                     edges = tibble::as_tibble(edges),
                                     directed = FALSE)




rivernetworkplotObj <- ggplot() +
  geom_sf(data = DRCprov, color = "#f0f0f0", fill = "#d9d9d9") +
  geom_sf(data = rivernetwork %>% activate(edges) %>% as_tibble() %>% st_as_sf(), color = "#1f78b4") +
  geom_sf(data = rivernetwork %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), size = 0.25, color = "#b2df8a") +
  geom_sf(data = ge, color = "#33a02c") +
  map_theme +
  theme(legend.position = "none") +
  coord_sf(datum = NA)

jpeg("results/figures/rivernetwork_clusters.jpg",
     height = 12, width = 8, units = "in", res = 300)
plot(rivernetworkplotObj)
graphics.off()


#............................................................................................................
# Get nearest neighbors for river network
#............................................................................................................
gecoords <- sf::st_coordinates(ge)

rivercoords <- rivernetwork %>%
  activate(nodes) %>%
  as_tibble() %>%
  st_as_sf() %>%
  sf::st_coordinates()

nn <- nabor::knn(data = rivercoords, query = gecoords, k=1)

dhsclust.dict <- ge$dhsclust
names(dhsclust.dict) <- unlist(nn$nn.idx)

# make tibble for search
dhsclust.tofrom <- tibble::as_tibble(t( combn(x = nn$nn.idx, m = 2) ))
colnames(dhsclust.tofrom) <- c("to", "from")



#............................................................................................................
# Get Length of Each Edge for Short Distance
# and Calculate Shortest Distance
#............................................................................................................
rivernetwork <- rivernetwork %>%
  tidygraph::activate(edges) %>%
  dplyr::mutate(length = sf::st_length(geometry))


#............................................................................................................
# Write out and run get shorest distance on LL
#............................................................................................................
saveRDS(rivernetwork, "data/derived_data/rivernetwork.rds")
saveRDS(dhsclust.tofrom, "data/derived_data/dhsclust.tofrom.rds")




# find dhsclusts from dictionary




find_replace_dictionary <- function(dict, replacevect){

  for(i in 1:nrow(dict)){
    replacevect[replacevect = names(dict)[i]] <- dict[i]
  }


}

saveRDS(dhsclust.tofrom,
        file = "data/distance_data/river_distance_forclusters.rds")





