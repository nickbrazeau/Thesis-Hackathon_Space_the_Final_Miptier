####################################################################################
## Purpose: Script to pull together covariates and genetic data
##
## Notes:
####################################################################################
#............................................................................
##### RASTER - Purpose of this script is to pull together covariate rasters ####
#...........................................................................
library(tidyverse)
library(sf)
library(raster)
source("R/basics.R")
source("R/pairwise_helpers.R")
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

# note assuming simple WGS84 here
#.............
# Pf Incidence
#.............
pfincidence <- raster::raster("data/raw_data/MAPrasters/getRaster/2019_Global_Pf_Incidence.201310_.18_40_8_2020_12_04.tiff")
pfincidence <- raster::mask(pfincidence, DRC)
pfincidence <- raster::aggregate(pfincidence, fact = 1.5, fun = mean)
summary(values(pfincidence))

#................................
# Urbanicity
#................................
urban <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_accessibility_to_cities_v1.0_latest_10_.18_40_8_2020_12_04.tiff")
urban <- raster::mask(urban, DRC)
values(urban)[values(urban) < -9998] <- NA
urban <- raster::aggregate(urban, fact = 6, fun = mean)
summary(values(urban))

#..........................................................
# Join
#..........................................................
# make raster stack
rstrscovars <- list(pfincidence, urban)


#..........................................................
# Let's make same resolution and extent
# want 0.05 which is urban, so will use that one
#..........................................................
rstrscovars.res <- lapply(rstrscovars, function(x){
  ret <- raster::projectRaster(x, urban)
  return(ret)
})

#..........................................................
# Save out Raw Format
#..........................................................
covar.rasterstack.raw <- raster::stack(rstrscovars.res)
names(covar.rasterstack.raw) <- c("incidence", "urban")
saveRDS(covar.rasterstack.raw, "data/derived_data/covar_rasterstack_raw.RDS")


#..........................................................
# Transformation of Values
#..........................................................
covar.rasterstack.derived <- covar.rasterstack.raw
# proportion back to real line
values(covar.rasterstack.derived$urban)[!is.na(values(covar.rasterstack.derived$urban))] <-  log(values(covar.rasterstack.derived$urban)[!is.na(values(covar.rasterstack.derived$urban))] + 1e-3)
values(covar.rasterstack.derived$incidence)[!is.na(values(covar.rasterstack.derived$incidence))] <-  logit(values(covar.rasterstack.derived$incidence)[!is.na(values(covar.rasterstack.derived$incidence))], tol = 1e-3)

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
#..............................................................
mtdt$urban_rura_ext[mtdt$hv001 %in% c("301")] <- 6e3
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




#................................................................................................................................
####### DISTANCES - Purpose of this section is to clean up distance measures and combine them #####
#................................................................................................................................
#.............................
# Cluster Level
#.............................
gcclust <- readRDS("data/distance_data/greater_circle_distance_forclusters.rds") %>%
  dplyr::mutate(item1 = as.character(item1),
                item2 = as.character(item2))
# roads
roadsclust <- readRDS("data/distance_data/clstr_road_distmeters_long.rds") %>%
  dplyr::mutate(item1 = as.character(item1),
                item2 = as.character(item2))

# airplanes
airclust <- readRDS("data/distance_data/clstr_airport_distmeters_long.rds") %>%
  dplyr::rename(item1 = hv001.x,
                item2 = hv001.y) %>%
  dplyr::mutate(item1 = as.character(item1),
                item2 = as.character(item2))

distancesmatrix_cluster <- dplyr::inner_join(gcclust, roadsclust, by = c("item1", "item2")) %>%
  dplyr::inner_join(., airclust, by = c("item1", "item2")) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y", "gcdistance", "roaddistance", "airportdistance"))


#................
# Out
#................
saveRDS(object = distancesmatrix_cluster, file = "data/distance_data/distancematrix_bycluster.rds")


#................................................................................................................................
#### GENETICS - Purpose of this section is to combine the genetic and mtdt data ####
#................................................................................................................................
#..............................................................
# Read in genetic dat
#..............................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.rds")
ibS <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS.long.rds")

#..............................................................
# Read in Metadata
#..............................................................
smplmtdt <- readRDS(file = "data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("name", "barcode", "hv001", "adm1name", "latnum", "longnum"))


#..............................................................
# Bring together ibD
#..............................................................
colnames(smplmtdt)[1] <- "smpl1"
ibD <- dplyr::left_join(ibD, smplmtdt, by = "smpl1")
colnames(smplmtdt)[1] <- "smpl2"
ibD <- dplyr::left_join(ibD, smplmtdt, by = "smpl2")

#..............................................................
# Bring together ibS
#..............................................................
colnames(smplmtdt)[1] <- "smpl1"
ibS <- dplyr::left_join(ibS, smplmtdt, by = "smpl1")
colnames(smplmtdt)[1] <- "smpl2"
ibS <- dplyr::left_join(ibS, smplmtdt, by = "smpl2")

#..............................................................
# Save out
#..............................................................
saveRDS(ibD,
        "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
saveRDS(ibS,
        "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS.long.mtdt.rds")




