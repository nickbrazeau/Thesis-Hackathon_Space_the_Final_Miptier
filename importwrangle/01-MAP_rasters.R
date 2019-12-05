#................................................................................................................................
# Purpose of this script is to download rasters from MAP
#................................................................................................................................

library(tidyverse)
library(malariaAtlas)
library(raster)
source("R/basics.R")
# bounding box
DRC <- sf::as_Spatial(osmdata::getbb("Democratic Republic of the Congo",
                                     featuretype = "country",
                                     format_out = 'sf_polygon'))
#................................
# Download Rasters
#................................

available_rasters <- malariaAtlas::listRaster()
rasterlist <- c("2019_Nature_Africa_Housing_2015", "2015_Nature_Africa_Incidence",
                "2015_Nature_Africa_PR", "2015_Nature_Africa_ITN", "2015_Nature_Africa_IRS",
                "2015_Nature_Africa_ACT", "2015_friction_surface_v1_Decompressed",
                "2015_accessibility_to_cities_v1.0")
raster.titles <- available_rasters %>%
  dplyr::filter(raster_code %in% rasterlist) %>%
  dplyr::filter(raster_code != "2019_Nature_Africa_Housing_2015") %>%
  dplyr::select("title") %>%
  unlist(.)

dir.create("data/raw_data/MAPrasters/", recursive = T)
lapply(raster.titles, function(x, bb){
  ret <- malariaAtlas::getRaster(surface = x, extent = bb,
                                 file_path = "data/raw_data/MAPrasters/")
  return(ret)}, bb = caf)

# have to do this one seperately (using a debug trick manually)
malariaAtlas::getRaster(surface = "Prevalence of improved housing",
                        extent = caf,
                        file_path = "data/raw_data/MAPrasters/")




#................................
# Read in and wrangle rasters
#................................
# download bigger to avoid NAs and then crop down
dir.create("data/derived_data/MAPrasters/")

# pr
parasiterate <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_Nature_Africa_PR_latest_10_.18_40_8_2019_11_14.tiff")%>%
  raster::mask(., DRC)
raster::writeRaster(x = parasiterate,
                    filename = "data/derived_data/MAPrasters/parasiterate.grd",
                    overwrite = T)


# act use
actuse <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_Nature_Africa_ACT_latest_10_.18_40_8_2019_11_14.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = actuse,
                    filename = "data/derived_data/MAPrasters/actuse.grd",
                    overwrite = T)


# net use
netuse <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_Nature_Africa_ITN_latest_10_.18_40_8_2019_11_14.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = netuse,
                    filename = "data/derived_data/MAPrasters/netuse.grd",
                    overwrite = T)

# house type
houseuse <- raster::raster("data/raw_data/MAPrasters/getRaster/2019_Nature_Africa_Housing_2015_latest_10_.18_40_8_2019_11_14.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = houseuse,
                    filename = "data/derived_data/MAPrasters/houseuse.grd",
                    overwrite = T)


# friction surface
frict <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_friction_surface_v1_Decompressed_latest_10_.18_40_8_2019_11_14.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = frict,
                    filename = "data/derived_data/MAPrasters/frict.grd",
                    overwrite = T)

# travel time
trav <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_accessibility_to_cities_v1.0_latest_10_.18_40_8_2019_11_14.tiff") %>%
  raster::mask(., DRC)
values(trav)[values(trav) == -9999] <- NA
raster::writeRaster(x = trav,
                    filename = "data/derived_data/MAPrasters/trav_acc.grd",
                    overwrite = T)




