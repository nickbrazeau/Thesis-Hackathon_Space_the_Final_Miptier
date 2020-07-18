#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import
# the world pop data from DRC 2013
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
# need this for bounding box
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

worldpop <- raster::raster("data/raw_data/worldpop/cod_ppp_2013.tif")
worldpop.drc <- raster::crop(x = worldpop, y = DRC)

# write out raw
dir.create("data/derived_data/worldpop/", recursive = T)
raster::writeRaster(worldpop.drc,  filename = "data/derived_data/worldpop/drc_worldpop_raw.grd",
                    overwrite=T)
