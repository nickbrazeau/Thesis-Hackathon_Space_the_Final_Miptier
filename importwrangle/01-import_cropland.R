#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle
# the night light composite data from MCD12Q1 2013
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sp)
library(sf)
source("R/basics.R")

# create bounding box for DRC
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

#..............................................................................
# Read in Land Coverage
#..............................................................................
landcov2013 <- raster::raster("data/raw_data/land_coverage/dataset-satellite-land-cover-902f410c-4e52-4fae-badd-2fc35ec86c62/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2013-v2.0.7cds.nc")
landcov2013.drc <- raster::crop(x = landcov2013, y = DRC)
#landcov2013.drc <- raster::projectRaster(from = landcov2013.drc, to = landcov2013.drc,
#                                         crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m

#..............................................................................
# Lift Over to Binary Cropland yes or no
#..............................................................................
# using data dictionary from here: https://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-QuickUserGuide-LC-Maps_v2-0-7.pdf
origvals <- raster::values(landcov2013.drc)
newvals <- ifelse(origvals %in% c(10, 20, 30, 40), 1, 0)
raster::values(landcov2013.drc) <- newvals

# write out
dir.create("data/derived_data/crops/", recursive = T)
raster::writeRaster(landcov2013.drc, filename = "data/derived_data/crops/cropland_surface.grd",
                    overwrite = T)
