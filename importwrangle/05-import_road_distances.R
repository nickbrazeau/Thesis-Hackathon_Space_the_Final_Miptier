library(tidyverse)
library(sf)
#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to  use OSRM to calculate road network distances
#     (1) All clusters
#     (2) All samples
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


diag(clstrtrvldists$distances) <- 1e6 # set this to really high to avoid self compare
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
  gc.err[i] <- geosphere::distHaversine(p1 = sf::as_Spatial(row.err[i,]),
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

#..................
# save out
#..................
clstrtrvldists.long <- clstrtrvldists.long %>%
  dplyr::rename(roaddistance = distance)
saveRDS(clstrtrvldists.long, file = "data/distance_data/clstr_road_distmeters_long.rds")





#.......................................................................................................
# Smpl Level
#.......................................................................................................
drcsmpls <- sf::st_as_sf(readRDS("data/distance_data/drcsmpls_foruse.rds"))
rownames(drcsmpls) <- drcsmpls$id
drcsmpls <- drcsmpls %>%
  dplyr::select(-c("hv001"))



# now access API
smplvldists <- osrm::osrmTable(loc = drcsmpls,
                               measure = "distance") # in meters


diag(smplvldists$distances) <- 1e6 # set this to really high to avoid self compare

#-----------------------------------------------------------------
# Process OSRM Outs
#-----------------------------------------------------------------
# same issues as above but duplicated since multiple obs per clust
which(smplvldists$distances == 0)
# note here that because we have duplicate clusters, some of these 0s are real
# but GC will also return 0 in that case, so better to be conservative and waste computational time
# figure out problem col and row
cols <- paste0("col", seq(1:ncol(smplvldists$distances)))
rows <- paste0("row", seq(1:nrow(smplvldists$distances)))
rowscols <- expand.grid(rows, cols)
rowscols <- paste0(rowscols[,1], "-", rowscols[,2])
dim(rowscols) <- dim(smplvldists$distances)

# misbehaving iterations
err <- rowscols[ which(smplvldists$distances == 0)]
err.df <- data.frame(dhsrow = stringr::str_split_fixed(err, "-", n=2)[,1],
                     dhscol = stringr::str_split_fixed(err, "-", n=2)[,2]) %>%
  dplyr::mutate(dhsrow = gsub("row", "", dhsrow),
                dhsrow = as.numeric(dhsrow),
                dhscol = gsub("col", "", dhscol),
                dhscol = as.numeric(dhscol)
  )

# Find responsible samples
err.df$dhsrow <- rownames(smplvldists$distances)[err.df$dhsrow]
err.df$dhscol <- colnames(smplvldists$distances)[err.df$dhscol]

# find the clusters that are responsible
row.err <- tibble::tibble(id = err.df$dhsrow) %>%
  dplyr::left_join(., y=drcsmpls, by = "id")
row.err <- sf::st_as_sf(row.err)

col.err <- tibble::tibble(id = err.df$dhscol) %>%
  dplyr::left_join(., y=drcsmpls, by = "id")
col.err <- sf::st_as_sf(col.err)


# Get GC Distance
gc.err <- rep(NA, nrow(row.err))
for(i in 1:nrow(row.err)){
  gc.err[i] <- geosphere::distHaversine(p1 = sf::as_Spatial(row.err[i,]),
                                        p2 = sf::as_Spatial(col.err[i,]))
}


# finally fix zeroes
smplvldists$distances[ which(smplvldists$distances == 0) ] <- gc.err


#..........
# Note, cluter 469 (with one sample) cannot be resolved with osrm
# do have one sample from that cluster...
# however cluster 313 is nearby, approximately 19763.25
# just going to add this approximate greater circle distance (2e4)
# to all 313 road distances to get 469 road dists
#..........
findreplacements.df <- sf::st_as_sf(readRDS("data/distance_data/drcsmpls_foruse.rds"))
sf::st_geometry(findreplacements.df) <- NULL
ge.replace <- ge %>%
  dplyr::rename(hv001 = dhsclust)
findreplacements.df <- dplyr::left_join(findreplacements.df, ge.replace, by = "hv001")

drcsmpls313 <- findreplacements.df %>%
  dplyr::filter(hv001 == 313)

drcsmpls469 <- findreplacements.df %>%
  dplyr::filter(hv001 == 469)

smplvldists.long <- broom::tidy(as.dist( smplvldists$distances ))

clstr313 <- smplvldists.long %>%
  dplyr::filter(item1 == "X2P7W" | item2 == "X2P7W") %>%
  dplyr::mutate(newitem1 = ifelse(item1 == "X2P7W", item1, item2),
                newitem2 = ifelse(item2 == "X2P7W", item1, item2)) %>%
  dplyr::arrange(newitem2)

clstr313$distance <- clstr313$distance + 2e4

# don't rearrange here
clstr469 <- smplvldists.long %>%
  dplyr::filter(item1 == "I3O3D" | item2 == "I3O3D") %>%
  dplyr::mutate(newitem2 = ifelse(item1 == "I3O3D", item2, item1) )


clstr469 <- left_join(clstr469, clstr313, by = "newitem2")

#..................
# now put dist back
#..................
smplvldists.long$distance[is.na(smplvldists.long$distance)] <- clstr469$distance.y

#..................
# save out
#..................
smplvldists.long <- smplvldists.long %>%
  dplyr::rename(roaddistance = distance)
saveRDS(smplvldists.long, file = "data/distance_data/smpl_road_distmeters_long.rds")
















