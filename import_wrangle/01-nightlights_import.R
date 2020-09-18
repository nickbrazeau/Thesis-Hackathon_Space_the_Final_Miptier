#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import
# the night light composite data from VIIRS
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sp)
library(sf)
source("R/basics.R")
# create bounding box of Central Africa for space
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"


#..............................................................................
# Night Light Raster Merge
#..............................................................................
# night light rasters
N0E60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_00N060E_v10_c201701311200/SVDNB_npp_20150101-20151231_00N060E_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N0W60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_00N060W_v10_c201701311200/SVDNB_npp_20150101-20151231_00N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N0W180 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_00N180W_v10_c201701311200/SVDNB_npp_20150101-20151231_00N180W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N75E60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_75N060E_v10_c201701311200/SVDNB_npp_20150101-20151231_75N060E_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N75W60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_75N060W_v10_c201701311200/SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N75W180 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_75N180W_v10_c201701311200/SVDNB_npp_20150101-20151231_75N180W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")

nightlights.rstrs <- list(N0E60, N0W60, N0W180,
                          N75E60, N75W60, N75W180)

nightlights.merge <- Reduce(function(...) merge(...), nightlights.rstrs)
# crop for size
nightlights.drc <- raster::crop(x = nightlights.merge, y = caf)

# check values
sum(values(nightlights.drc) < 0)
# apparently these <0 values can result from too much correction https://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=5888
# given that there are so few, I am going to set them to NA
values(nightlights.drc)[values(nightlights.drc) < 0] <- NA

# write out raw
dir.create("data/derived_data/nightlights/", recursive = T)
raster::writeRaster(nightlights.drc,  filename = "data/derived_data/nightlights/drc_nightlights_raw.grd",
                    overwrite=T)


