#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data from HOTOSM (open street map)
# primarily interested health sites and waterways
#----------------------------------------------------------------------------------------------------
gdrive <- file.choose()

# libraries
library(tidyverse)
library(sf)
library(rgeos)
tol = 1e-3
# read in GE as import
ge <- sf::st_as_sf(readRDS(paste0(gdrive, "/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

# Note manually downloading these from site
# https://data.humdata.org/dataset/

#----------------------------------------------------------------------------------------------------
# Waterways
#----------------------------------------------------------------------------------------------------
wtrlns <- sf::read_sf(paste0(gdrive, "/data/raw_data/hotosm_data/hotosm_cod_waterways_lines_shp/hotosm_cod_waterways_lines.shp")) %>%
  dplyr::filter(waterway %in% c("stream", "river", "riverbank")) %>%
  dplyr::rename(water = waterway) %>%
  dplyr::select(c("osm_id", "water", "geometry"))

wtrply <- sf::read_sf(paste0(gdrive, "/data/raw_data/hotosm_data/hotosm_cod_waterways_polygons_shp/hotosm_cod_waterways_polygons.shp")) %>%
  dplyr::filter(water == "lake") %>%
  dplyr::select(c("osm_id", "water", "geometry"))

wtr <- sf::st_combine(rbind(wtrlns, wtrply))
wtr <-  sf::st_union( wtr )

wtrdist <- sf::st_distance(x = ge,
                           y = wtr)
wtrdist_out <- data.frame(
  hv001 = ge$dhsclust,
  wtrdist_cont_clst = apply(wtrdist, 1, min)
)
#----------------------------------------------------------------------------------------------------
# Health Sites
#----------------------------------------------------------------------------------------------------
hlthsites <- sf::read_sf(paste0(gdrive, "/data/raw_data/hotosm_data/hotosm_drc_healthsites_shapefiles/healthsites.shp")) %>%
  dplyr::filter(type %in% c("hospital", "clinic"))
htlhdist <- sf::st_distance(x = ge,
                            y = hlthsites)

hlthdist_out <- data.frame(
  hv001 = ge$dhsclust,
  hlthdist_cont_clst = apply(htlhdist, 1, min)
)



#---------------------------------------------------------------------------------
# Precipation Data from CHIRPS
#---------------------------------------------------------------------------------
heavyRain::getCHIRPS(region = "africa",
                     format = "tifs",
                     tres = "monthly",
                     sres = 0.05, # near same resolution as Manny pulled down
                     begin = as.Date("2013-01-01"),
                     end = as.Date("2014-12-31"),
                     dsn = paste0(gdrive, "/data/raw_data/weather_data/CHIRPS/"),
                     overwrite = T)




#----------------------------------------------------------------------------------------------------
# write out
#----------------------------------------------------------------------------------------------------
system('gunzip $gdrive/data/raw_data/weather_data/CHIRPS/*') # need a system loc here
saveRDS(object = wtrdist_out, file = paste0(gdrive, "/data/derived_data/hotosm_waterways_dist.rds"))
saveRDS(object = hlthdist_out, file = paste0(gdrive, "/data/derived_data/hotosm_healthsites_dist.rds"))
