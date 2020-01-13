library(tidyverse)
library(sf)
library(raster)
#---------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to  use OSRM to calculate road network distances
#     (1) for clusters
#     (2) All provs
#---------------------------------------------------------------------------------------------------------------------------------------
# drc prov for plot
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))

#.......................................................................................................
# Cluster Level
#.......................................................................................................
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::select(c("hv001"))

#..................................
# SPIN up Docker and start server
# then add these options
#..................................
# https://github.com/rCarto/osrm
remotes::install_github("rCarto/osrm")
library(osrm)
options(osrm.server = "http://0.0.0.0:5000/", osrm.profile = "driving")

ge.osrm <- ge %>%
  dplyr::select("hv001")
rownames(ge.osrm) <- ge.osrm$hv001

#..............................................................
# pull in friction surface
#..............................................................
# This is what the surface actually is
# https://explorer.earthengine.google.com/#detail/Oxford%2FMAP%2Ffriction_surface_2015_v1_0
frict <- raster::raster("data/derived_data/MAPrasters/frict.grd")
frict <- raster::projectRaster(from = frict, to = frict,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m"))

jpeg("results/figures/MAP_friction_surface.jpg", width = 8, height = 11, res = 400, units = "in")
ggplot() +
  ggspatial::layer_spatial(data = frict, aes(fill = stat(band1))) +
  prettybasemap_nodrc_dark +
  geom_sf(data = DRCprov, color = "#737373", fill = NA, size = 0.05) +
  ggtilte("MAP Minutes to travel One Meter")
graphics.off()

#..............................................................
# Make cluster by cluster comparison grid
#..............................................................
clst <- ge$hv001[!duplicated(ge$hv001)]
clst.comb <- as.data.frame(t(combn(clst, 2))) %>%
  magrittr::set_colnames(c("source", "destination"))


router <- function(source, destination){
  #..............................................................
  # Find Routes
  #..............................................................
  src <- ge %>% dplyr::filter(hv001 == source)
  dst <- ge %>% dplyr::filter(hv001 == destination)

  ret <- osrm::osrmRoute(src = src, dst = dst,
                         overview = "full", returnclass = "sf")
  distance <- osrm::osrmTable(src = src, dst = dst,
                              measure = "distance")
  distance <- as.vector(distance$distances)

  if (is.null(ret)) {
    return(NA)
  } else {

    cells <- unlist( raster::extract(x = frict,
                           y = ret) )
    # 2.512473
    # about 1km by 1 km
  }
  x = unique(cells)
  w <- as.vector(table(cells))
  Wimean <- stats::weighted.mean(x = x, w = w)

  ret <- list(
    distance = distance,
    duration = distance * Wimean
  )


  return(ret)
}

clst.comb$routes <- purrr::pmap(clst.comb, router)

nick = purrr::pmap(clst.comb[which(clst.comb$source == 469), ], router)

#-----------------------------------------------------------------
# Process OSRM Outs
#-----------------------------------------------------------------

#..........
# Note, cluter 469 cannot be resolved with osrm
# do have one sample from that cluster...
# however cluster 313 is nearby, approximately 19763.25 meters away
#..........



#..............................................................
# Extract Cells now and sum
# Note, we are in meters/minute
#..............................................................

duration <- raster::extract()









