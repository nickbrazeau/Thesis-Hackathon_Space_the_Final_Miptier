#................................................................................................................................
# Purpose of this script is to calculate cluster duration
# in minutes to the nearest hospital
#................................................................................................................................
library(tidyverse)
library(sf)


# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

hlthsites.harvard.drc <- readxl::read_excel("data/raw_data/harvard_dataverse/Ouma_Okiro_Snow_Africa_Hospitals_Data.xlsx") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  magrittr::set_colnames(gsub(pattern = " ", "_", colnames(.))) %>%
  dplyr::filter(country == "Democratic Republic of Congo") %>%
  sf::st_as_sf(coords = c("long", "lat"),
               crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))


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


hlthsites.harvard.drc.osrm <- hlthsites.harvard.drc %>%
  dplyr::mutate(id = paste0("hlth", seq(1:nrow(.))),
                id = factor(id)) %>%
  dplyr::select(id)

rownames(hlthsites.harvard.drc.osrm) <- hlthsites.harvard.drc.osrm$id


# now access API
hosptrvltmes <- osrm::osrmTable(src = ge.osrm,
                                dst = hlthsites.harvard.drc.osrm,
                                measure = "duration") # in minutes


#-----------------------------------------------------------------
# Process OSRM Outs
#-----------------------------------------------------------------

#..........
# Note, 4 Clusters have 0s that indicate that OSRM could
# not resolve a route between the two points via a road
# For an example, spin up the OSRM frontend and look at these two points:
# S1.2523, E19.1959  --- S1.1085, E19.0815
# Basically what happens is that we don't consider "non-road" paths. So both these points
# touch a small segment of the road and that's it. Similar to a triangle but we don't count
# the distances of the lines which is non-road (so the distance is 0)
# These zero cells are going to be replaced by
# greater circle distances * (50km/hr * 60min/hr) as an "off-roading" pace
#..........
# r counts matrices by columns
which(hosptrvltmes$durations == 0)
# figure out problem col and row
cols <- paste0("col", seq(1:ncol(hosptrvltmes$durations)))
rows <- paste0("row", seq(1:nrow(hosptrvltmes$durations)))
rowscols <- expand.grid(rows, cols)
rowscols <- paste0(rowscols[,1], "-", rowscols[,2])
dim(rowscols) <- dim(hosptrvltmes$durations)

# misbehaving iterations
err <- rowscols[ which(hosptrvltmes$durations == 0)]
err.df <- data.frame(dhsrow = stringr::str_split_fixed(err, "-", n=2)[,1],
                     hlthcol = stringr::str_split_fixed(err, "-", n=2)[,2]) %>%
  dplyr::mutate(dhsrow = gsub("row", "", dhsrow),
                dhsrow = as.numeric(dhsrow),
                hlthcol = gsub("col", "", hlthcol),
                hlthcol = as.numeric(hlthcol)
  )

# Find responsible points
err.df$dhsrow <- rownames(hosptrvltmes$durations)[err.df$dhsrow]
# note, hlthzones are just "hlth" paste0 with rownum (see above), dhsclust are not rownums because of lost clusters
# extract points, note repeat clusters
ge.err <- tibble::tibble(dhsclust = as.numeric(err.df$dhsrow)) %>%
  dplyr::left_join(., y=ge, by = "dhsclust")
ge.err <- sf::st_as_sf(ge.err)

hlthsites.err <- hlthsites.harvard.drc[err.df$hlthcol, ]
# sanity
nrow(ge.err) == nrow(hlthsites.err)
# Get GC Distance
gc.err <- rep(NA, nrow(ge.err))
for(i in 1:nrow(ge.err)){
  gc.err[i] <- geosphere::distHaversine(p1 = sf::as_Spatial(ge.err[i,]),
                                        p2 = sf::as_Spatial(hlthsites.err[i,]))
}

# now convert meters to minutes
gc.err <- gc.err * (60/(1e3*50)) # min per meter (assuming we travel off-road at 40kph)

# finally fix zeroes
hosptrvltmes$durations[ which(hosptrvltmes$durations == 0) ] <- gc.err



hlthdist_out <-
  cbind.data.frame(hv001 = as.numeric( rownames(hosptrvltmes$duration) ), hosptrvltmes$duration ) %>%
  tidyr::gather(., key = "hsptl", value = "hlthst_duration", 2:ncol(.)) %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    hlthst_nrst_duration = min(hlthst_duration)
  )


#..........
# Note, cluter 469 cannot be resolved with osrm
# however cluster 313 is nearby, approximately 19763.25
# just going to add this greater circle distance
#..........
clst469 <- ge %>%
  dplyr::filter(dhsclust == 469)

dist <- raster::pointDistance(p1 = sf::as_Spatial(clst469),
                              p2 = sf::as_Spatial(ge),
                              lonlat = T)
# find nearest cluster
nghbrfor469 <- ge$dhsclust[which(dist %in%  sort(dist)[2])]
hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == nghbrfor469]
sort(dist)[2] # need to travel about 20 km, which is about 1hr by car in rural setting

hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == 469] <- hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == nghbrfor469] + 60




#..............
# out
#..............
saveRDS(object = hlthdist_out, file = "data/derived_data/hlthdist_out_minduration.rds")
