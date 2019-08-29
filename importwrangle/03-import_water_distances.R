#................................................................................................................................
# Purpose of this script is to import that waterways
# need to calculate distance by water
#................................................................................................................................

# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))


#----------------------------------------------------------------------------------------------------
# Waterways
#----------------------------------------------------------------------------------------------------
wtrlns <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_lines_shp/hotosm_cod_waterways_lines.shp") %>%
  dplyr::filter(waterway %in% c("stream", "river", "riverbank")) %>%
  dplyr::rename(watertype = waterway) %>%
  dplyr::select(c("osm_id", "watertype", "geometry"))

wtrply <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_polygons_shp/hotosm_cod_waterways_polygons.shp") %>%
  dplyr::filter(water == "lake") %>%
  dplyr::rename(watertype = water) %>%
  dplyr::select(c("osm_id", "watertype", "geometry"))

wtr <- sf::st_combine(rbind(wtrlns, wtrply))
wtr <-  sf::st_union( wtr )

wtrdist <- sf::st_distance(x = ge,
                           y = wtr)
wtrdist_out <- data.frame(
  hv001 = ge$dhsclust,
  wtrdist_cont_clst = apply(wtrdist, 1, min)
)


save(wtrply, wtrlns, file = "data/derived_data/hotosm_waterways.rda")











#----------------------------------------------------------------------------------------------------
# Do something with distance
#----------------------------------------------------------------------------------------------------


pts <- as.matrix(init[,c("longnum", "latnum")])[,1:2] # odd behavior bc it keep geometry ... not like a dataframe

t <- shp2graph::points2network(ntdata = sf::as_Spatial(road),
                               pointsxy = points,
                               approach = 1) #https://rdrr.io/cran/shp2graph/man/points2network.html
t.cnt <- nt.connect(sf::as_Spatial(road))
plot(t.cnt)
u <- shp2graph::nel2igraph(nodelist = t[[1]], edgelist = t[[2]])
plot(u)

igraph::shortest_paths(u, from = t)
shortest.paths(u, v=V(u), to=V(u))

v <- tidygraph::as_tbl_graph(u)
