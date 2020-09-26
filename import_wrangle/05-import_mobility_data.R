#########################################################################
# Purpose:
#
# Author: Nicholas F. Brazeau
#
# Data Source: World Pop Internal Flows https://www.worldpop.org/geodata/summary?id=1281
#       Methods: https://www.nature.com/articles/sdata201666
#                https://academic.oup.com/migration/article/3/1/89/2413406
#########################################################################
library(tidyverse)
library(sf)
source("R/basics.R")
# assume simple planar geo WGS84

# create bounding box of Central Africa for envelop
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
caf <- sf::st_bbox(sf::st_as_sf(caf))
caf <- matrix(
  c(caf["xmin"], caf["ymin"],
    caf["xmin"], caf["ymax"],
    caf["xmax"], caf["ymax"],
    caf["xmax"], caf["ymin"],
    caf["xmin"], caf["ymin"]),
  ncol = 2, byrow = T
)
caf <- sf::st_polygon(list(caf))

# drc polygon for intersection
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))

#............................................................
# read in data
#...........................................................
mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("name", "barcode", "hv001", "longnum", "latnum"))

flows <- readr::read_csv("data/raw_data/worldpop/COD_5yrs_InternalMigFlows_2010/COD_5yrs_InternalMigFlows_2010.csv")
centroids <- sf::read_sf("data/raw_data/worldpop/COD_5yrs_InternalMigFlows_2010/COD_AdminUnit_Centroids/COD_AdminUnit_Centroids.shp")
# ggplot checek centroids
ggplot() +
  geom_sf(data = centroids, color = "red") +
  ggrepel::geom_text_repel(data = centroids, aes(x = POINT_X, y = POINT_Y, label = IPUMSID))
# some very closely overlapping points...

#......................
# voronoi tesselations from migration data
#......................
# https://stackoverflow.com/questions/45719790/create-voronoi-polygon-with-simple-feature-in-r
vr <- sf::st_voronoi(sf::st_union(centroids), envelope = caf)

# quick sanity
ggplot() +
  geom_sf(data = vr) +
  geom_sf(data = DRCprov, alpha = 0.25) +
  geom_sf(data = centroids, color = "red")

# clip to DRC polygons
vr <- st_intersection(st_cast(vr), st_union(DRCprov))
# sanity check
plot(st_intersection(st_cast(vr), st_union(DRCprov)), col = 0)
plot(centroids, add = TRUE)



#............................................................
# bring together
#...........................................................
# take back to dataframe
vrdf <- sf::st_as_sf(vr) %>%
  dplyr::mutate(ID = 1:nrow(.))

#......................
# get nearest neighbor from original
#......................
newcentroids <- vrdf %>%
  dplyr::mutate(centroid = sf::st_centroid(x))

ipusmid <- split(centroids, 1:nrow(centroids))
ipusmid_match_list <- lapply(ipusmid, function(x){gcdist <- sf::st_distance(x, newcentroids$centroid)
matchid <- newcentroids$ID[which(gcdist == min(gcdist))]
ret <- tibble::tibble(IPUMSID = x$IPUMSID,
                      ID = matchid)
return(ret)
}) %>%
  dplyr::bind_rows()

length(unique(ipusmid_match_list$IPUMSID))
length(unique(ipusmid_match_list$ID))
#......................
# for three sites, where original points were close together (20, 21, 25) Knn not as clear
#......................
newcentroidsplotdf <- newcentroids %>%
  dplyr::mutate(long = sf::st_coordinates(centroid)[,1],
                lat = sf::st_coordinates(centroid)[,2]) %>%
  dplyr::select(-c("x", "centroid"))
st_geometry(newcentroidsplotdf) <- NULL

origcentroidsplotdf <- centroids  %>%
  dplyr::mutate(POINT_X = sf::st_coordinates(geometry)[,1],
                POINT_Y = sf::st_coordinates(geometry)[,2]) %>%
  dplyr::select(-c("geometry"))
st_geometry(origcentroidsplotdf) <- NULL


plotdf <- dplyr::left_join(vrdf, ipusmid_match_list, by = "ID") %>%
  dplyr::left_join(., newcentroidsplotdf, by = "ID") %>%
  dplyr::left_join(., origcentroidsplotdf, by = "IPUMSID")

ggplot() +
  geom_sf(data = vrdf) +
  geom_point(data = plotdf, aes(x = long, y = lat), color = "red") +
  geom_point(data = plotdf, aes(x = POINT_X, y = POINT_Y), color = "blue") +
  ggrepel::geom_text_repel(data = plotdf, aes(x = long, y = lat, label = ID), color = "red") +
  ggrepel::geom_text_repel(data = plotdf, aes(x = POINT_X, y = POINT_Y, label = IPUMSID), color = "blue")


#......................
# manually adjust some duplicates
#......................
dups <- ipusmid_match_list$ID[duplicated(ipusmid_match_list$ID)]
goodmatches <- ipusmid_match_list %>%
  dplyr::filter(! ID %in% dups)
badmatches <- ipusmid_match_list %>%
  dplyr::filter(ID %in% dups)
badmatches$ID[badmatches$IPUMSID == 1] <- 21
badmatches$ID[badmatches$IPUMSID == 38] <- 20
badmatches$ID[badmatches$IPUMSID == 5] <- 25
# bring together
ipusmid_match_list <- dplyr::bind_rows(goodmatches, badmatches)


#......................
# bring it together
#......................
vrdf <- vrdf %>%
  dplyr::left_join(., ipusmid_match_list, by = "ID") %>%
  dplyr::select(-c("ID"))
saveRDS(vrdf, "data/distance_data/voroni_base.RDS")


#............................................................
# Connect with Flows
#...........................................................
vrdf_from <- vrdf %>%
  dplyr::rename(NODEI = IPUMSID,
                poly_from = x)
vrdf_to <- vrdf %>%
  dplyr::rename(NODEJ = IPUMSID,
                poly_to = x)

flows <- flows %>%
  dplyr::left_join(., vrdf_from, by = "NODEI")
flows <- flows %>%
  dplyr::left_join(., vrdf_to, by = "NODEJ")

flows <- flows %>%
  dplyr::mutate(PrdMIG_scaled = my.scale(PrdMIG))

#......................
# bring it together
#......................
saveRDS(flows, "data/distance_data/voroni_migration_flows_fromworldpop.RDS")

