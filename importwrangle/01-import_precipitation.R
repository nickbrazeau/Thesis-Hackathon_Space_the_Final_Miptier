#---------------------------------------------------------------------------------
# Purpose of this script is to download precipation data
# and find the lagged precipitation per cluster
#---------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(raster)
source("R/basics.R")
readRasterBB <- function(rstfile, bb = bb){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = bb)
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m

  return(ret)

}


#.............................................................
# Download Precipitation (CHRIPS)
#.............................................................
dir.create("data/raw_data/weather_data/CHIRPS/", recursive = T)
heavyRain::getCHIRPS(region = "africa",
                     format = "tifs",
                     tres = "monthly",
                     sres = 0.05,
                     begin = as.Date("2013-01-01"),
                     end = as.Date("2014-12-31"),
                     dsn = "data/raw_data/weather_data/CHIRPS/",
                     overwrite = T)

system('gunzip data/raw_data/weather_data/CHIRPS/*')





#.............................................................
# Date Wrangle
#.............................................................
dt <- readRDS("data/derived_data/DHS_qPCR_allkids_geo.rds")
dt <- dt %>%
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")))

# NOTE, some clusters have survey start and end dates that are in two months
# (eg boundaries aren't clean/coinciding with a month. Naturally). Given
# grouping by month, need to assign a clusters "month" on the majority of days
# that were spent surveying that clusters

# clusters without clean boundaries
clst_mnth_bounds <- dt[, c("hv001", "hvdate_dtdmy")] %>%
  dplyr::mutate(mnth = lubridate::month(hvdate_dtdmy)) %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(moremnths = length(unique(mnth))) %>%
  dplyr::filter(moremnths > 1)

clst_mnth_bounds.assign <- dt[, c("hv001", "hvdate_dtdmy")] %>%
  dplyr::filter(hv001 %in% clst_mnth_bounds$hv001) %>%
  dplyr::mutate(hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")) %>%
  dplyr::group_by(hv001, hvyrmnth_dtmnth) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == max(n)) %>%
  dplyr::select(-c("n"))


dates <- readr::read_csv("internal_datamap_files/pr_date_liftover.csv")
dt <- dt %>%
  dplyr::left_join(x=., y = clst_mnth_bounds.assign, by = "hv001") %>%
  dplyr::mutate(hvyrmnth_dtmnth = ifelse(is.na(hvyrmnth_dtmnth),
                                         paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-"),
                                         hvyrmnth_dtmnth))

dates <- readr::read_csv("internal_datamap_files/pr_date_liftover.csv")
dt <- dt %>%
  dplyr::left_join(x=., y=dates, by = "hvyrmnth_dtmnth") %>%
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))

xtabs(~dt$hvyrmnth_dtmnth + dt$hvyrmnth_dtmnth_lag)
# subset dt
dt <- dt %>%
  dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag", "urban_rura", "geometry"))

#.............................................................
# Manipulate Precipitation (CHRIPS)
#.............................................................

# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"

# precip data
precip <- list.files(path = "data/raw_data/weather_data/CHIRPS/", full.names = T,
                     pattern = ".tif")
precipfrst <- lapply(precip, readRasterBB, bb = caf)

precipdf <- tibble::tibble(names = basename(precip)) %>%
  dplyr::mutate(names = gsub("chirps-v2.0.|.tif", "", names),
                year = stringr::str_split_fixed(names, "\\.", n=2)[,1] ,
                mnth =  stringr::str_split_fixed(names, "\\.", n=2)[,2] ,
                hvdate_dtdmy = lubridate::dmy(paste(1, mnth, year, sep = "/")),
                year = lubridate::year(hvdate_dtdmy),
                mnth = lubridate::month(hvdate_dtdmy),
                hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
                precipraster = precipfrst) %>%
  dplyr::select(c("hvyrmnth_dtmnth_lag", "precipraster"))



wthrnd <- dt[,c("hv001", "hvyrmnth_dtmnth_lag", "geometry", "urban_rura")] %>%
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
wthrnd <- wthrnd[!duplicated(wthrnd$hv001),]

wthrnd <- wthrnd %>%
  dplyr::left_join(., precipdf)

# Drop in a for loop
wthrnd$precip_lag_cont_clst <- NA

for(i in 1:nrow(wthrnd)){
  # precip
  wthrnd$precip_lag_cont_clst[i] <-
    raster::extract(x = wthrnd$precipraster[[i]],
                    y = sf::as_Spatial(wthrnd$geometry[i]),
                    buffer = wthrnd$buffer[i],
                    fun = mean,
                    sp = F
    )


}

wthrnd <- wthrnd %>%
  dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag", "precip_lag_cont_clst")) %>%
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))
sf::st_geometry(wthrnd) <- NULL
wthrnd <- wthrnd %>%
  dplyr::mutate(
    precip_lag_cont_scale_clst = my.scale(precip_lag_cont_clst, center = T, scale = T)
  )


saveRDS(object = wthrnd, file = "data/derived_data/CHIRPS_cluster_precip.rds")

