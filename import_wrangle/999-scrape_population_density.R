#................................................................................................................................
# Purpose of this script is to download and
# scrape population density from worldpop
#................................................................................................................................
library(tidyverse)
library(raster)
source("R/basics.R")

# create bounding box for DRC mask
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
# create bounding box of Central Africa for download speed
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"

#............................................................
# DHS Approach
#...........................................................
dhs <- readRDS("data/raw_data/dhsdata/datasets/CDGC62FL.rds") %>%
  dplyr::rename(hv001 = DHSCLUST) %>%
  dplyr::select(c("hv001", "All_Population_Count_2015"))

# space mips clusters
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
clsts <- unique(c(ibD$hv001.x, ibD$hv001.y))

# down select
dhs <- dhs %>%
  dplyr::filter(hv001 %in% clsts) %>%
  dplyr::rename(dhs_pop_dens = All_Population_Count_2015)
summary(dhs$dhs_pop_dens)
hist(dhs$dhs_pop_dens)

#............................................................
# Raster Approach
#...........................................................
#................................
# Read in Rasters
#................................
worldpop <- raster::raster("data/raw_data/worldpop/cod_ppp_2013.tif")

#................................
# Extract out for density at sampling locations
#................................
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

popdensdf <- sf::st_as_sf(mtdt, coords = c("longnum", "latnum"), crs = 4326)

#............................................................
# scrape out population densities
#...........................................................
#popdensdf$urban_rura_ext[popdensdf$hv001 %in% c("301")] <- 6e3
popdensdf$popdens <- NA

for (i in 1:nrow(popdensdf)){
  popdensdf$popdens[i] <- raster::extract(x = worldpop,
                                       y = sf::as_Spatial(popdensdf$geometry[i]),
                                       buffer = popdensdf$urban_rura_ext[i],
                                       fun = sum,
                                       na.rm = T)
}
#......................
# look at output
#......................
worldpop.agg <- raster::aggregate(worldpop, fact = 100, fun = sum)
summary(popdensdf$popdens)
hist(popdensdf$popdens)
popdensdf %>%
  dplyr::left_join(., mtdt) %>%
  ggplot() +
  ggspatial::layer_spatial(data = worldpop.agg,
                           aes(fill = stat(band1))) +
    geom_point(aes(x = longnum, y = latnum, color = popdens)) +
  scale_fill_viridis_c()

popdensdf %>%
  dplyr::left_join(., mtdt) %>%
  dplyr::filter(popdens == 0) %>%
  ggplot() +
  ggspatial::layer_spatial(data = worldpop.agg,
                           aes(fill = stat(band1))) +
  geom_point(aes(x = longnum, y = latnum, color = popdens)) +
  scale_fill_viridis_c()

#......................
# save out
#......................
popdensdf <- popdensdf %>%
  dplyr::select(c("hv001", "popdens")) %>%
  dplyr::rename(worldpop_pop_dens = popdens) %>%
  dplyr::left_join(., dhs, by = "hv001")
saveRDS(popdensdf, "data/derived_data/popdens_samplinglocations_raw.RDS")




