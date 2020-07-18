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
rasterlist <- c("2015_Nature_Africa_ITN", "2015_Nature_Africa_IRS",
                "2015_Nature_Africa_ACT", "2015_friction_surface_v1_Decompressed",
                "2015_accessibility_to_cities_v1.0")
raster.titles <- available_rasters %>%
  dplyr::filter(raster_code %in% rasterlist) %>%
  #dplyr::filter(raster_code != "2019_Nature_Africa_Housing_2015") %>%
  dplyr::select("title") %>%
  unlist(.)

dir.create("data/raw_data/MAPrasters/", recursive = T)
raster.years <- rep(NA, length(raster.titles))
raster.years[grepl("2019_Global_Pf_Incidence", raster.titles)] <- 2013

malariaAtlas::getRaster(surface = raster.titles,
                        year = raster.years,
                        shp = caf,
                        file_path = "data/raw_data/MAPrasters/")


# have to do this one seperately (just failed above)
malariaAtlas::getRaster(surface = "Prevalence of improved housing 2015",
                        shp = caf,
                        file_path = "data/raw_data/MAPrasters/")

malariaAtlas::getRaster(surface = "Plasmodium  falciparum  Incidence.",
                        year = 2013,
                        shp = caf,
                        file_path = "data/raw_data/MAPrasters/")

# malariaAtlas::getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015",
#                         shp = caf,
#                         file_path = "data/raw_data/MAPrasters/")


#................................
# Read in and wrangle rasters
#................................
# download bigger to avoid NAs and then crop down
dir.create("data/derived_data/MAPrasters/")

# pr
pfincidence <- raster::raster("data/raw_data/MAPrasters/getRaster/2019_Global_Pf_Incidence.201310_.18_40_8_2020_03_31.tiff")%>%
  raster::mask(., DRC)
raster::writeRaster(x = pfincidence,
                    filename = "data/derived_data/MAPrasters/pfincidence.grd",
                    overwrite = T)
plot(pfincidence)

# act use
actuse <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_Nature_Africa_ACT_latest_10_.18_40_8_2020_03_31.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = actuse,
                    filename = "data/derived_data/MAPrasters/actuse.grd",
                    overwrite = T)
plot(actuse) # this raster is useless for DRC

# net use
netuse <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_Nature_Africa_ITN_latest_10_.18_40_8_2020_03_31.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = netuse,
                    filename = "data/derived_data/MAPrasters/netuse.grd",
                    overwrite = T)
plot(netuse)

# house type
houseuse <- raster::raster("data/raw_data/MAPrasters/getRaster/2019_Nature_Africa_Housing_2015_latest_10_.18_40_8_2020_03_31.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = houseuse,
                    filename = "data/derived_data/MAPrasters/houseuse.grd",
                    overwrite = T)
plot(houseuse)

# friction surface
frict <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_friction_surface_v1_Decompressed_latest_10_.18_40_8_2020_03_31.tiff") %>%
  raster::mask(., DRC)
raster::writeRaster(x = frict,
                    filename = "data/derived_data/MAPrasters/frict.grd",
                    overwrite = T)
plot(frict)

# travel time
trav <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_accessibility_to_cities_v1.0_latest_10_.18_40_8_2020_03_31.tiff") %>%
  raster::mask(., DRC)
values(trav)[values(trav) == -9999] <- NA
raster::writeRaster(x = trav,
                    filename = "data/derived_data/MAPrasters/trav_acc.grd",
                    overwrite = T)
plot(trav)



