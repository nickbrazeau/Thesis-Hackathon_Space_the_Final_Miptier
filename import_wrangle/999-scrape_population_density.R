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

#................................
# Read in Rasters
#................................
worldpop <- raster::raster("data/raw_data/worldpop/cod_ppp_2013.tif")


#..............................................................
# Extract out for density at sampling locations
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

popdensdf <- sf::st_as_sf(mtdt, coords = c("longnum", "latnum"), crs = 4326)

#............................................................
# scrape out population densities
#...........................................................
#popdensdf$urban_rura_ext[popdensdf$hv001 %in% c("301")] <- 6e3
popdensdf$popdens <- NA

for (i in 1:nrow(popdensdf)){
  popdensdf$popdens <- raster::extract(x = worldpop,
                                       y = sf::as_Spatial(popdensdf$geometry[i]),
                                       buffer = popdensdf$urban_rura_ext[i],
                                       fun = sum,
                                       na.rm = T)
}
#......................
# look at output
#......................
summary(popdensdf$popdens)
hist(popdensdf$popdens)
popdensdf %>%
  ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(aes(x = longnum, y = latnum, color = popdens))
#......................
# save out
#......................
popdensdf <- popdensdf %>%
  dplyr::select(c("hv001", "popdens"))
saveRDS(popdensdf, "data/derived_data/popdens_samplinglocations_raw.RDS")




