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
# Waterways and simplify
#...............................
drcrivers <- sf::st_read("data/raw_data/river_data/combined/combind_rivers_postgrass/combined_rivers_postgrass.shp")
drcrivers.simp <- shp2graph::nt.connect(sf::as_Spatial(drcrivers))
drcrivers.simp <- sf::st_as_sf(drcrivers.simp)

# need to add back in geometry
sf::st_crs(drcrivers.simp) <- sf::st_crs(drcrivers)


# save this out for plotting later
saveRDS(object = drcrivers.simp,
        file = "data/derived_data/river_network.RDS")


#............................................................................................................
# Manipulate Shapes to Prepare for Network
#............................................................................................................
#...............................
# Edges
# Give each line a unique ID
#...............................
edges <- drcrivers.simp %>%
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

edges <- edges %>%
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
  prettybasemap_nodrc_dark +
  geom_sf(data = DRCprov, fill = "#525252", color = "#737373") +
  geom_sf(data = rivernetwork %>% activate(edges) %>% as_tibble() %>% st_as_sf(),
          color = "#9ecae1", size = 0.05) +
  geom_sf(data = rivernetwork %>% activate(nodes) %>% as_tibble() %>% st_as_sf(),
          size = 0.05, color = "#9ecae1") +
  geom_sf(data = ge, color = "#ff2e2e") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']),
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']),
           datum = NA)

jpeg("results/figures/rivernetwork_clusters.jpg",
     height = 12, width = 8, units = "in", res = 300)
plot(rivernetworkplotObj)
graphics.off()
saveRDS(rivernetworkplotObj, file = "data/distance_data/rivernetworkplotObj.RDS")

#............................................................................................................
# Get nearest neighbors for river network
#............................................................................................................
gecoords <- ge %>%
  dplyr::mutate(
    X = sf::st_coordinates(ge)[,1],
    Y = sf::st_coordinates(ge)[,2]
    ) %>%
  dplyr::select(c("dhsclust", "X", "Y"))
sf::st_geometry(gecoords) <- NULL

rivercoords <- rivernetwork %>%
  activate(nodes) %>%
  as_tibble() %>%
  st_as_sf() %>%
  sf::st_coordinates()

nn <- nabor::knn(data = rivercoords, query = gecoords[,c("X", "Y")], k=1)


# make tibble for search
dhsclust.tofrom <- tibble::as_tibble(t( combn(x = gecoords$dhsclust, m = 2) ))
rivernet.tofrom <- tibble::as_tibble(t( combn(x = nn$nn.idx, m = 2) ))
dhsclust.tofrom <- cbind.data.frame(dhsclust.tofrom, rivernet.tofrom)
colnames(dhsclust.tofrom) <- c("hv001.x", "hv001.y", "to", "from")



#............................................................................................................
# Get Length of Each Edge for Short Distance
# and Calculate Shortest Distance
#............................................................................................................
rivernetwork <- rivernetwork %>%
  tidygraph::activate(edges) %>%
  dplyr::mutate(length = sf::st_length(geometry))


#............................................................................................................
# get shorest distance
#............................................................................................................
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

dhsclust.tofrom$riverdist <- furrr::future_pmap(dhsclust.tofrom[,c("to", "from")],
                                                get_shortest_distance_length)
dhsclust.tofrom.unnested <- dhsclust.tofrom %>%
  tidyr::unnest(cols = riverdist)


#............................................................................................................
# Liftover to get cluster hv001 labels
#............................................................................................................


saveRDS(dhsclust.tofrom,
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
# Get nearest neighbors for Province and run shortest distance calculation
#............................................................................................................
drcpovcoords <- drcpov %>%
  dplyr::mutate(
    X = sf::st_coordinates(drcpov)[,1],
    Y = sf::st_coordinates(drcpov)[,2]
    ) %>%
  dplyr::select(c("adm1name", "X", "Y"))
sf::st_geometry(drcpovcoords) <- NULL

nn <- nabor::knn(data = rivercoords, query = drcpovcoords[,c("X", "Y")],  k=1)



# make tibble for search
drcpovcoords.tofrom <- tibble::as_tibble(t( combn(x = drcpovcoords$adm1name, m = 2) ))
rivernetprov.tofrom <- tibble::as_tibble(t( combn(x = nn$nn.idx, m = 2) ))
drcpovcoords.tofrom <- cbind.data.frame(drcpovcoords.tofrom, rivernetprov.tofrom)
colnames(drcpovcoords.tofrom) <- c("adm1name.x", "adm1name.y", "to", "from")

# run function
drcpovcoords.tofrom$riverdist <- purrr::pmap(drcpovcoords.tofrom[,c("to", "from")],
                                             get_shortest_distance_length)

drcpovcoords.tofrom <- drcpovcoords.tofrom %>%
  tidyr::unnest(cols = riverdist)

#..............................................................
# LIFTOVER to from table to prov names
#..............................................................
drcpovcoords.tofrom <- drcpovcoords.tofrom %>%
  dplyr::select(c("adm1name.x", "adm1name.y", "riverdist")) %>%
  dplyr::rename(from = adm1name.x,
                to = adm1name.y)


saveRDS(drcpovcoords.tofrom,
        file = "data/distance_data/river_distance_forprovinces.rds")


