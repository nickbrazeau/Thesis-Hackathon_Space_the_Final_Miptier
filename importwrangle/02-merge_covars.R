#................................................................................................................................
# Purpose of this script is to pull together covariate rasters
#................................................................................................................................
library(raster)
source("R/basics.R")
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
#.............
# Pf Prevalence
#.............
parasiterate <- raster::raster("data/derived_data/MAPrasters/parasiterate.grd")

#................................
# Malaria Interventions
#................................
actuse <- raster::raster("data/derived_data/MAPrasters/actuse.grd")
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
  parasiterate, # prev
  precip, temp, elev, # ecological
  crops, actuse, netuse, houseuse, # malaria interv/risk
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
names(covar.rasterstack.raw) <- c("prev", "precip", "temp", "elev",
                                  "crops", "actuse", "netuse", "housing",
                                  "urban")

# force crops to binary
values(covar.rasterstack.raw$crops)[!values(covar.rasterstack.raw$crops) %in% c(0,1)] <- round(values(covar.rasterstack.raw$crops)[!values(covar.rasterstack.raw$crops) %in% c(0,1)], 0)
saveRDS(covar.rasterstack.raw, "data/derived_data/covar_rasterstack_raw.RDS")


#..........................................................
# Transformation of Values
#..........................................................
covar.rasterstack.derived <- covar.rasterstack.raw
# proportion back to real line
values(covar.rasterstack.derived$prev) <-  logit(values(covar.rasterstack.derived$prev), tol = 1e-3)

# weather and elev
values(covar.rasterstack.derived$precip) <- my.scale(values(covar.rasterstack.derived$precip))
values(covar.rasterstack.derived$temp) <- my.scale(values(covar.rasterstack.derived$temp))
values(covar.rasterstack.derived$elev) <- my.scale(values(covar.rasterstack.derived$elev))


# crops are binary, no transform

# proportion back to real line
values(covar.rasterstack.derived$actuse) <- my.scale( logit(values(covar.rasterstack.derived$actuse), tol = 1e-3) )
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


drcsmpls <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/drcsmpls_foruse.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("id", "hv001")) %>%
  dplyr::rename(name = id)

mtdt <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(name = id) %>%
  dplyr::filter(name %in% drcsmpls$name) %>%
  dplyr::select(c("hv001", "longnum", "latnum")) %>%
  dplyr::left_join(., ge) %>%
  dplyr::filter(!duplicated(.))

mtdt.sf <- sf::st_as_sf(mtdt, coords = c("longnum", "latnum"), crs = 4326)



#.........................
# Raw Covariates
#.........................
covar.raw.extract <- matrix(NA, nrow = nrow(mtdt.sf),
                            ncol = length(names(covar.rasterstack.raw)))
for (i in 1:nrow(mtdt.sf)){
  covar.raw.extract[i, ] <- raster::extract(x = covar.rasterstack.raw,
                                            y = sf::as_Spatial(mtdt.sf$geometry[i]),
                                            buffer = mtdt$urban_rura_ext[i],
                                            fun = mean,
                                            na.rm = T)
}

# 301 clusters has missing values around it -- impute mean
impvalues <- apply(covar.raw.extract, 2, function(x){
  if (any(is.na(x))) {
    ret <- mean(x, na.rm = T)
    return(ret)
  }
})

covar.raw.extract[is.na(covar.raw.extract)] <- unlist(impvalues)



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


# 301 clusters has missing values around it -- impute mean
impvalues <- apply(covar.scaled.extract, 2, function(x){
  if (any(is.na(x))) {
    ret <- mean(x, na.rm = T)
    return(ret)
  }
})

covar.scaled.extract[is.na(covar.scaled.extract)] <- unlist(impvalues)


covar.scaled.extract <- cbind.data.frame(hv001 = mtdt.sf$hv001, covar.scaled.extract)
colnames(covar.scaled.extract) <- c("hv001", names(covar.rasterstack.derived))
saveRDS(covar.scaled.extract, "data/derived_data/covar_rasterstack_samplinglocations_scaled.RDS")

#..............................................................
# Extract out for covariates for prov
#..............................................................
DRCprov <- sf::st_as_sf(readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/map_bases/gadm/gadm36_COD_1_sp.rds"))

#.........................
# Raw Covariates
#.........................
# set up items for extraction
covar.raw.extract <- matrix(NA, nrow = nrow(DRCprov),
                            ncol = length(names(covar.rasterstack.raw)))

# fill in table
for (i in 1:nrow(DRCprov)){
  covar.raw.extract[i, ] <-raster::extract(x = covar.rasterstack.raw,
                                           y = sf::as_Spatial(DRCprov$geometry[i]),
                                           fun = mean,
                                           na.rm = T,
                                           sp = F)}

covar.raw.extract <- cbind.data.frame(adm1name = DRCprov$adm1name,
                                      covar.raw.extract)
colnames(covar.raw.extract) <- c("adm1name", names(covar.rasterstack.raw))
saveRDS(covar.raw.extract, "data/derived_data/covar_rasterstack_provlocations_raw.RDS")

