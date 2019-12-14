library(tidyverse)
library(sf)
#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to  use OSRM to calculate road network distances
#     (1) for clusters
#     (2) All provs
#---------------------------------------------------------------------------------------------------------------------------------------

#.......................................................................................................
# Cluster Level
#.......................................................................................................
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


diag(clstrtrvldists$distances) <- .Machine$double.xmax # set this to really high to avoid self compare
#-----------------------------------------------------------------
# Process OSRM Outs
#-----------------------------------------------------------------

#..........
# Note, 2 pairs of clusters have 0s that indicate that OSRM could
# not resolve a route between the two points via a road
# For an example, spin up the OSRM frontend and look at these two points:
# S1.2523, E19.1959  --- S1.1085, E19.0815
# Basically what happens is that we don't consider "non-road" paths. So both these points
# touch a small segment of the road and that's it. Similar to a triangle but we don't count
# the distances of the lines which is non-road (so the distance is 0)
# These zero cells are going to be replaced by
# greater circle distances
#..........
# r counts matrices by columns
which(clstrtrvldists$distances == 0)
# figure out problem col and row
cols <- paste0("col", seq(1:ncol(clstrtrvldists$distances)))
rows <- paste0("row", seq(1:nrow(clstrtrvldists$distances)))
rowscols <- expand.grid(rows, cols)
rowscols <- paste0(rowscols[,1], "-", rowscols[,2])
dim(rowscols) <- dim(clstrtrvldists$distances)

# misbehaving iterations
err <- rowscols[ which(clstrtrvldists$distances == 0)]
err.df <- data.frame(dhsrow = stringr::str_split_fixed(err, "-", n=2)[,1],
                     dhscol = stringr::str_split_fixed(err, "-", n=2)[,2]) %>%
  dplyr::mutate(dhsrow = gsub("row", "", dhsrow),
                dhsrow = as.numeric(dhsrow),
                dhscol = gsub("col", "", dhscol),
                dhscol = as.numeric(dhscol)
  )

# Find responsible points
err.df$dhsrow <- rownames(clstrtrvldists$distances)[err.df$dhsrow]
err.df$dhscol <- colnames(clstrtrvldists$distances)[err.df$dhscol]


row.err <- tibble::tibble(dhsclust = as.numeric(err.df$dhsrow)) %>%
  dplyr::left_join(., y=ge, by = "dhsclust")
row.err <- sf::st_as_sf(row.err)

col.err <- tibble::tibble(dhsclust = as.numeric(err.df$dhscol)) %>%
  dplyr::left_join(., y=ge, by = "dhsclust")
col.err <- sf::st_as_sf(col.err)


# Get GC Distance
gc.err <- rep(NA, nrow(row.err))
for(i in 1:nrow(row.err)){
  gc.err[i] <- geosphere::distGeo(p1 = sf::as_Spatial(row.err[i,]),
                                  p2 = sf::as_Spatial(col.err[i,]))
}


# finally fix zeroes
clstrtrvldists$distances[ which(clstrtrvldists$distances == 0) ] <- gc.err


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

# finally, find 313 and 469 comparison and set to 2e4
clstrtrvldists.long$distance[clstrtrvldists.long$item1 == 469 & clstrtrvldists.long$item2 == 313] <- 2e4


#..................
# save out
#..................
clstrtrvldists.long <- clstrtrvldists.long %>%
  dplyr::rename(roaddistance = distance)

saveRDS(clstrtrvldists.long, file = "data/distance_data/clstr_road_distmeters_long.rds")



#.......................................................................................................
# Province Level
#.......................................................................................................
# read in DRCprov as import
drcpov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")) %>%
  dplyr::mutate(provcentroid = sf::st_centroid(geometry))


drcpov.osrm <- drcpov %>%
  dplyr::select(c("adm1name", "provcentroid"))
# make centroid geom
sf::st_geometry(drcpov.osrm) <- NULL
drcpov.osrm <- sf::st_as_sf(drcpov.osrm) # now point is geom


# altering the centroid for the Mai-Ndombe Province, which has an odd shape is
# having a hard time being resolved by OSRM
# to -2.529588, 19.050577
crs <- sf::st_crs(drcpov.osrm$provcentroid[drcpov.osrm$adm1name == "Mai-Ndombe"])
drcpov.osrm$provcentroid[drcpov.osrm$adm1name == "Mai-Ndombe"] <- sf::st_sfc(st_point(c(19.3506,
                                                                                        -2.5240)),
                                                                   crs = crs)



# now access API
provlvldists <- osrm::osrmTable(loc = drcpov.osrm,
                                  measure = "distance") # in meters



#..................
# save out
#..................
provlvldists.long <- provlvldists %>%
  .$distances %>%
  as.dist(.) %>%
  broom::tidy(.) %>%
  dplyr::rename(roaddistance = distance)

saveRDS(provlvldists.long, file = "data/distance_data/prov_road_distmeters_long.rds")


#..............................................................
# Make plot Obj for Road Network
#..............................................................
# README: I exported out the pdf lines portion of the
# congo-democratic-republic-latest lines and to a shape file
# and am reading it in below
roads <- sf::st_read("data/raw_data/osrm/congo-democratic-republic-latest_shape_LINES_fromPbf.shp")

roads.filt <- roads %>%
  dplyr::filter( highway %in% c(
    "primary", "secondary",
    "primary_link", "secondary_link",
    "tertiary", "tertiary_link",
    "trunk", "trunk_link", "motorway", "road"
  ))

#...............................
# map import
#...............................
load("~/Documents/GitHub/Space_the_Final_Miptier/data/map_bases/space_mips_maps_bases.rda")
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

# plot
roadnetworkplotObj <- ggplot() +
  prettybasemap_nodrc_dark +
  geom_sf(data = DRCprov, fill = "#525252", color = "#737373", size = 0.05) +
  geom_sf(data = roads.filt, fill = "#f0f0f0", color = "#f0f0f0") +
  geom_sf(data = ge, color = "#ff2e2e") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']),
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']),
           datum = NA)

jpeg("results/figures/roadnetwork_clusters.jpg",
     height = 12, width = 8, units = "in", res = 300)
plot(roadnetworkplotObj)
graphics.off()
saveRDS(roadnetworkplotObj, file = "data/distance_data/roadnetworkplotObj.RDS")

