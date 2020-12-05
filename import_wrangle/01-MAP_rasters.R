#................................................................................................................................
# Purpose of this script is to download rasters from MAP
#................................................................................................................................
library(tidyverse)
remotes::install_github("malaria-atlas-project/malariaAtlas")
library(malariaAtlas)
library(raster)
source("R/basics.R")

# create bounding box for DRC mask
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
# create bounding box of Central Africa for download speed
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"

#................................
# Download Rasters
#................................
available_rasters <- malariaAtlas::listRaster()
raster.titles <- available_rasters %>%
  dplyr::filter(raster_code %in% "2015_accessibility_to_cities_v1.0") %>%
  dplyr::select("title") %>%
  unlist(.)

dir.create("data/raw_data/MAPrasters/", recursive = T)
malariaAtlas::getRaster(surface = raster.titles,
                        shp = caf,
                        file_path = "data/raw_data/MAPrasters/")

malariaAtlas::getRaster(surface = "Plasmodium  falciparum  Incidence.",
                        year = 2013,
                        shp = caf,
                        file_path = "data/raw_data/MAPrasters/") # specify year





