#................................................................................................................................
# Purpose of this script is to calculate cluster duration
# in minutes to the nearest hospital
#................................................................................................................................
library(tidyverse)
library(sf)


# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

hlthsites.harvard.drc <- readxl::read_excel("data/raw_data/harvard_dataverse/Ouma_Okiro_Snow_Africa_Hospitals_Data.xlsx") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  magrittr::set_colnames(gsub(pattern = " ", "_", colnames(.))) %>%
  dplyr::filter(country == "Democratic Republic of Congo") %>%
  sf::st_as_sf(coords = c("long", "lat"),
               crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))


#..................................
# SPIN up Docker and start server
# then add these options
#..................................
# https://github.com/rCarto/osrm
remotes::install_github("rCarto/osrm")
library(osrm)
options(osrm.server = "http://0.0.0.0:5000/", osrm.profile = "driving")

ge.osrm <- ge %>%
  dplyr::select("dhsclust")
rownames(ge.osrm) <- ge.osrm$dhsclust


hlthsites.harvard.drc.osrm <- hlthsites.harvard.drc %>%
  dplyr::mutate(id = paste0("hlth", seq(1:nrow(.))),
                id = factor(id)) %>%
  dplyr::select(id)

rownames(hlthsites.harvard.drc.osrm) <- hlthsites.harvard.drc.osrm$id


# now access API
hosptrvltmes <- osrm::osrmTable(src = ge.osrm,
                                dst = hlthsites.harvard.drc.osrm,
                                measure = "duration") # in minutes


hlthdist_out <-
  cbind.data.frame(hv001 = as.numeric( rownames(hosptrvltmes$duration) ), hosptrvltmes$duration ) %>%
  tidyr::gather(., key = "hsptl", value = "hlthst_duration", 2:ncol(.)) %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    hlthst_nrst_duration = min(hlthst_duration)
  )



#..........
# Note, cluter 469 cannot be resolved with osrm
# however cluster 313 is nearby, approximately 19763.25
# just going to add this greater circle distance
#..........
clst469 <- ge %>%
  dplyr::filter(dhsclust == 469)

dist <- raster::pointDistance(p1 = sf::as_Spatial(clst469),
                              p2 = sf::as_Spatial(ge),
                              lonlat = T)
# find nearest cluster
nghbrfor469 <- ge$dhsclust[which(dist %in%  sort(dist)[2])]
hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == nghbrfor469]
sort(dist)[2] # need to travel about 20 km, which is about 1hr by car in rural setting

hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == 469] <- hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == nghbrfor469] + 60

#..............
# out
#..............
saveRDS(object = hlthdist_out, file = "data/derived_data/hlthdist_out_minduration.rds")
