#---------------------------------------------------------------------------------
# Purpose of this script
# is to download an elevation raster for the DRC
#---------------------------------------------------------------------------------
library(raster)
library(elevatr)

# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

# raster to ggplot for color
dem.raster <- elevatr::get_elev_raster(DRC, z=4)  # elevatr expects sp, similar res to MAP

# assume that values less than 0 are not land in the DRC people inhabit
# (people aren't living below sea level)
values(dem.raster)[values(dem.raster) < 0 ] <- NA
# save out this surface
dir.create("data/derived_data/elevation/", recursive = T)
raster::writeRaster(dem.raster,  filename = "data/derived_data/elevation/drc_elevation.grd",
                    overwrite = T)


