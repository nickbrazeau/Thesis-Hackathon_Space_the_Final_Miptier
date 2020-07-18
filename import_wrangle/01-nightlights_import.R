#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle
# the night light composite data from VIIRS
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sp)
library(sf)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
# need this for bounding box
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

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
nightlights.drc <- raster::crop(x = nightlights.merge, y = DRC)
#nightlights.drc <- raster::projectRaster(from = nightlights.drc, to = nightlights.drc,
#                                         crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m


# check values
sum(values(nightlights.drc) < 0) # there are 2237 values less than 0 out of 44920800
# apparently these <0 values can result from too much correction https://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=5888
# given that there are so few, I am going to set them to NA
values(nightlights.drc)[values(nightlights.drc) < 0] <- NA

# write out raw
dir.create("data/derived_data/nightlights/", recursive = T)
raster::writeRaster(nightlights.drc,  filename = "data/derived_data/nightlights/drc_nightlights_raw.grd",
                    overwrite=T)


# look at data for clusters
summary(values(nightlights.drc))
hist( values(nightlights.drc) )
hist( values(nightlights.drc)[values(nightlights.drc) > 0] )


