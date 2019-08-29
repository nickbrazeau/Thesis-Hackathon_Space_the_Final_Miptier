#................................................................................................................................
# Purpose of this script is to calculate
# greater cirlce distance between clusters
#................................................................................................................................

ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))
ge.nosf <- ge
sf::st_geometry(ge.nosf) <- NULL


gc <- ge.nosf %>%
  dplyr::select(c("longnum", "latnum")) %>%
  geosphere::distm(x =., fun = geosphere::distGeo)



gc.long <- broom::tidy(as.dist( gc ))



saveRDS(gc.long, file = "data/distance_data/greater_circle_distance.rds")
