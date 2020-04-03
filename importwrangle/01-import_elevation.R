#---------------------------------------------------------------------------------
# Purpose of this script
# is to download an elevation raster for the DRC
#---------------------------------------------------------------------------------
library(raster)
library(elevatr)

# bounding box
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

# raster to ggplot for color
dem.raster <- elevatr::get_elev_raster(DRC, z=4)  # elevatr expects sp, similar res to MAP
dem.raster <- raster::mask(dem.raster, DRC)
# save out this surface
dir.create("data/derived_data/elevation/", recursive = T)
raster::writeRaster(dem.raster,  filename = "data/derived_data/elevation/drc_elevation.grd",
                    overwrite = T)


