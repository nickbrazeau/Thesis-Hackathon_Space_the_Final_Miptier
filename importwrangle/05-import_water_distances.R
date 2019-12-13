#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate the water distanc between
# (1) clusters
# (2) Provinces
# Following this tutorial https://www.r-spatial.org/r/2019/09/26/spatial-networks.html
# to get shorest distances along the river network
# Above allows us to make a tidy spatial network
#---------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(tidygraph)
library(igraph)
library(raster)
source("R/themes.R")
source("R/pairwise_helpers.R")

#............................................................................................................
# IMPORT DATA for CLUSTER SECTION
#............................................................................................................
#...............................
# Ge import
#...............................
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))



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



load("~/Documents/GitHub/Space_the_Final_Miptier/data/map_bases/space_mips_maps_bases.rda")
rivernetworkplotObj <- ggplot() +
  geom_sf(data = DRCprov, fill = "#525252", color = "#737373") +
  geom_sf(data = rivernetwork %>% activate(edges) %>% as_tibble() %>% st_as_sf(),
          color = "#9ecae1") +
  geom_sf(data = rivernetwork %>% activate(nodes) %>% as_tibble() %>% st_as_sf(),
          size = 0.25, color = "#9ecae1") +
  geom_sf(data = ge, color = "#ff2e2e") +
  prettybasemap_nodrc_dark +
  theme(legend.position = "none")

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
saveRDS(dhsclust.tofrom, "data/distance_data/dhsclust.tofrom.rds")

#............................................................................................................
# Read in from LL
#............................................................................................................
LLdhsclust.tofrom <- readRDS(file = "data/distance_data/indexed_river_dist_long.rds")

fromjoin = tibble::tibble(from = as.numeric( names(dhsclust.dict) ), dhsclustfrom = dhsclust.dict)
tojoin = tibble::tibble(to = as.numeric( names(dhsclust.dict) ), dhsclustto = dhsclust.dict)

LLdhsclust.tofrom <- LLdhsclust.tofrom %>%
  dplyr::left_join(., fromjoin, by = "from") %>%
  dplyr::left_join(., tojoin, by = "to")

LLdhsclust.tofrom <- LLdhsclust.tofrom[!duplicated(LLdhsclust.tofrom), ] # nearest neighbor result can be duplicate for multiple clusters


saveRDS(LLdhsclust.tofrom,
        file = "data/distance_data/river_distance_forclusters.rds")



#............................................................................................................
# IMPORT DATA for PROVINCE SECTION
#............................................................................................................
# read in DRCprov as import
drcpov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")) %>%
  dplyr::mutate(provcentroid = sf::st_centroid(geometry))


drcpov <- drcpov %>%
  dplyr::select(c("adm1name", "provcentroid"))
# make centroid geom
sf::st_geometry(drcpov) <- NULL
drcpov <- sf::st_as_sf(drcpov) # make centroids dominating geom

#............................................................................................................
# Get nearest neighbors for river network
#............................................................................................................
drcpovcoords <- sf::st_coordinates(drcpov)

nn <- nabor::knn(data = rivercoords,
                 query = drcpovcoords,
                 k=1)

drcpovcoords.dict <- drcpov$adm1name
names(drcpovcoords.dict) <- unlist(nn$nn.idx)

# make tibble for search
drcpovcoords.tofrom <- tibble::as_tibble(t( combn(x = nn$nn.idx, m = 2) ))
colnames(drcpovcoords.tofrom) <- c("to", "from")



#............................................................................................................
# Get Length of Each Edge for Short Distance
# and Calculate Shortest Distance
#............................................................................................................
rivernetwork <- rivernetwork %>%
  tidygraph::activate(edges) %>%
  dplyr::mutate(length = sf::st_length(geometry))



get_shortest_distance_length <- function(to, from){
  path <- igraph::shortest_paths(
    graph = rivernetwork,
    from = from,
    to = to,
    output = 'both',
    weights = rivernetwork %>% activate(edges) %>% pull(length))

  dist <- subgraph.edges(rivernetwork, eids = path$epath %>% unlist()) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    as_tibble() %>%
    summarise(length = sum(length)) %>%
    dplyr::pull(length)

  return(dist)
}

# run function
drcpovcoords.tofrom$riverdist <- purrr::pmap(drcpovcoords.tofrom,
                                             get_shortest_distance_length)

drcpovcoords.tofrom <- drcpovcoords.tofrom %>%
  tidyr::unnest(cols = riverdist)

#..............................................................
# LIFTOVER to from table to prov names
#..............................................................

fromjoin = tibble::tibble(from = as.numeric( names(drcpovcoords.dict) ),
                          dhsprovfrom = drcpovcoords.dict)
tojoin = tibble::tibble(to = as.numeric( names(drcpovcoords.dict) ),
                        dhsprovto = drcpovcoords.dict)

drcpovcoords.tofrom <- drcpovcoords.tofrom %>%
  dplyr::left_join(., fromjoin, by = "from") %>%
  dplyr::left_join(., tojoin, by = "to")

# liftover
# prov dict
provname.dict <- tibble::tibble(item1 = 1:26,
                                dhsprovfrom = drcpov$adm1name)
drcpovcoords.tofrom <- dplyr::left_join(drcpovcoords.tofrom, provname.dict,
                                        by = "dhsprovfrom")

provname.dict <- tibble::tibble(item2 = 1:26,
                                dhsprovto = drcpov$adm1name)
drcpovcoords.tofrom <- dplyr::left_join(drcpovcoords.tofrom, provname.dict,
                                        by = "dhsprovto")

drcpovcoords.tofrom <- drcpovcoords.tofrom[!duplicated(drcpovcoords.tofrom), ] # nearest neighbor result can be duplicate for multiple clusters


saveRDS(drcpovcoords.tofrom,
        file = "data/distance_data/river_distance_forprovinces.rds")


