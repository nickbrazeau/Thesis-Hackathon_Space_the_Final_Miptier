#............................................................................................................
# Purpose of this script is to run on cluster and get
# the nearest point along a line segment with respect to each dhsclust
#............................................................................................................
library(sf)
library(tidyverse)
library(tidygraph)
library(igraph)


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


#..................................................
# Read in Data and run main
#..................................................
rivernetwork <- readRDS("data/derived_data/rivernetwork.rds")
dhsclust.tofrom <- readRDS("data/distance_data/dhsclust.tofrom.rds")


dhsclust.tofrom$riverdist <- purrr::pmap(dhsclust.tofrom,
                                         get_shortest_distance_length)

dhsclust.tofrom <- dhsclust.tofrom %>%
  tidyr::unnest(cols = riverdist)

saveRDS(object = dhsclust.tofrom,
        file = "data/distance_data/indexed_river_dist_long.rds")



