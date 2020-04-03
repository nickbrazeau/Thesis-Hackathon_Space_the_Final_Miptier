#................................................................................................................................
# Purpose of this script is to pull together covariate rasters
#................................................................................................................................
library(tidyverse)
library(sf)
library(raster)
source("R/basics.R")
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

#.............
# Pf Incidence
#.............
pfincidence <- raster::raster("data/derived_data/MAPrasters/pfincidence.grd")

#................................
# Malaria Interventions
#................................
netuse <- raster::raster("data/derived_data/MAPrasters/netuse.grd")
houseuse <- raster::raster("data/derived_data/MAPrasters/houseuse.grd")

#................................
# Urbanicity
#................................
urban <- raster::raster("data/derived_data/urbanicity_raster/urbanicity.grd")

#.............
# Weather
#.............
precip <- raster::raster("data/derived_data/weather/precipitation.grd")
temp <- raster::raster("data/derived_data/weather/temperature.grd")

#.............
# Cluster-Level Altitude
#.............
elev <- raster::raster("data/derived_data/elevation/drc_elevation.grd")

#.............
# Cropland
#.............
crops <- raster::raster("data/derived_data/crops/cropland_surface.grd")
crops <- raster::aggregate(crops, fact = 18, fun = mean)


#..........................................................
# Join
#..........................................................
# make raster stack
rstrscovars <- list(
  pfincidence, # prev
  precip, temp, elev, # ecological
  crops, netuse, houseuse, # malaria interv/risk
  urban
)


#..........................................................
# Let's make same resolution and extent
# want 0.05 which is precip, so will use that one
#..........................................................
rstrscovars.res <- lapply(rstrscovars, function(x){
  ret <- raster::projectRaster(x,
                               precip)
  return(ret)
  })


#..........................................................
# Save out Raw Format
#..........................................................
covar.rasterstack.raw <- raster::stack(rstrscovars.res)
names(covar.rasterstack.raw) <- c("incidence", "precip", "temp", "elev",
                                  "crops", "netuse", "housing",
                                  "urban")

# force crops to binary
values(covar.rasterstack.raw$crops)[!values(covar.rasterstack.raw$crops) %in% c(0,1)] <- round(values(covar.rasterstack.raw$crops)[!values(covar.rasterstack.raw$crops) %in% c(0,1)], 0)
saveRDS(covar.rasterstack.raw, "data/derived_data/covar_rasterstack_raw.RDS")


#..........................................................
# Transformation of Values
#..........................................................
covar.rasterstack.derived <- covar.rasterstack.raw
# proportion back to real line
values(covar.rasterstack.derived$incidence) <-  logit(values(covar.rasterstack.derived$incidence), tol = 1e-3)

# weather and elev
values(covar.rasterstack.derived$precip) <- my.scale(values(covar.rasterstack.derived$precip))
values(covar.rasterstack.derived$temp) <- my.scale(values(covar.rasterstack.derived$temp))
values(covar.rasterstack.derived$elev) <- my.scale(values(covar.rasterstack.derived$elev))


# crops are binary, no transform

 # proportion back to real line
values(covar.rasterstack.derived$netuse) <- my.scale( logit(values(covar.rasterstack.derived$netuse), tol = 1e-3) )
# proportion back to real line
values(covar.rasterstack.derived$housing) <- my.scale( logit(values(covar.rasterstack.derived$housing), tol = 1e-3) )


# note urbanicity from PCA so already scaled

#..........................
#  save out
#..........................
saveRDS(covar.rasterstack.derived, "data/derived_data/covar_rasterstack_derived.RDS")

#..............................................................
# Extract out for covariates at sampling locations
#..............................................................
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::select(c("hv001", "urban_rura")) %>%
  dplyr::mutate(urban_rura_ext = ifelse(urban_rura == "R", 10000, 2000))
sf::st_geometry(ge) <- NULL


mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("hv001", "longnum", "latnum")) %>%
  dplyr::left_join(., ge, by = "hv001") %>%
  dplyr::filter(!duplicated(.))

mtdt.sf <- sf::st_as_sf(mtdt, coords = c("longnum", "latnum"), crs = 4326)


#.........................
# Raw Covariates
#.........................
#..............................................................
# NOTE, PREVIOUSLY
#  cluster 301 is missing incidence
#  cluster 301 and 223 are missing precip
#  cluster 223 is missing temp
#  cluster 301 is missing elevation
#  cluster 301 is missing netuse
#  cluster 301 is missing housing
#       both boundaries that likely have to do with cut and are urban
#       To fix, extend mean search to 10km from 2km
#..............................................................

mtdt$urban_rura_ext[mtdt$hv001 %in% c("223", "301")] <- 6e3

covar.raw.extract <- matrix(NA, nrow = nrow(mtdt.sf),
                            ncol = length(names(covar.rasterstack.raw)))
for (i in 1:nrow(mtdt.sf)){
  covar.raw.extract[i, ] <- raster::extract(x = covar.rasterstack.raw,
                                            y = sf::as_Spatial(mtdt.sf$geometry[i]),
                                            buffer = mtdt$urban_rura_ext[i],
                                            fun = mean,
                                            na.rm = T)
}

covar.raw.extract <- cbind.data.frame(hv001 = mtdt.sf$hv001, covar.raw.extract)
colnames(covar.raw.extract) <- c("hv001", names(covar.rasterstack.raw))
saveRDS(covar.raw.extract, "data/derived_data/covar_rasterstack_samplinglocations_raw.RDS")

#.........................
# Scaled Covariates
#.........................
covar.scaled.extract <- matrix(NA, nrow = nrow(mtdt.sf),
                            ncol = length(names(covar.rasterstack.derived)))
for (i in 1:nrow(mtdt.sf)){
  covar.scaled.extract[i, ] <- raster::extract(x = covar.rasterstack.derived,
                                            y = sf::as_Spatial(mtdt.sf$geometry[i]),
                                            buffer = mtdt.sf$urban_rura_ext[i],
                                            fun = mean,
                                            na.rm = T)
}

covar.scaled.extract <- cbind.data.frame(hv001 = mtdt.sf$hv001, covar.scaled.extract)
colnames(covar.scaled.extract) <- c("hv001", names(covar.rasterstack.derived))
saveRDS(covar.scaled.extract, "data/derived_data/covar_rasterstack_samplinglocations_scaled.RDS")


