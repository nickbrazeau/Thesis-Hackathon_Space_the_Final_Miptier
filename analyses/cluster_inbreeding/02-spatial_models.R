#########################################################################
# Purpose: Spatial Modeling of Cluster Inbreeding Coeffs
#
# Author: Nicholas F. Brazeau
#
# Date: April 08 2020
#########################################################################
source("R/basics.R")
source("R/gauss_proc_simple_functions.R")
library(tidyverse)
library(PrevMap)

#..............................................................
# read in data
#..............................................................
clst_inbd <- readRDS("data/derived_data/clst_inbreeding_dat/clst_inbreeding_convfilt_results.RDS")
clst_inbd <- clst_inbd %>%
  tidyr::gather(., key = "hv001", value = "Fclst", 3:ncol(.))
clst_inbd_m <- clst_inbd %>%
  dplyr::filter(hv001 == "m") %>%
  dplyr::rename(param = hv001)
clst_inbd_f <- clst_inbd %>%
  dplyr::filter(hv001 != "m" & hv001 != "param_set")

#......................
# split F clust results
#......................
clst_inbd_f <- clst_inbd_f %>%
  dplyr::group_by(param_set, spacetype) %>%
  tidyr::nest()

#......................
# geo
#......................
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))
clsts <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("hv001", "latnum", "longnum")) %>%
  dplyr::mutate(hv001 = as.character(hv001)) %>%
  dplyr::filter(!duplicated(.))



temp <- dplyr::left_join(clst_inbd_f$data[[1]], clsts, by = "hv001")
temp <- temp %>%
  dplyr::mutate(Fclst_logit = logit(as.numeric(Fclst))) %>%
  dplyr::filter(!hv001 %in% c("m", "conv_pass"))

#..............................................................
# temp
#..............................................................

poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14))
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("long","lat")

Fclst_raster <- fit_pred_spMLE(data = temp,
                               outcome = "Fclst_logit", covar = "1",
                               long_var = "longnum", lat_var = "latnum",
                               grid.pred = grid.pred, kappa = 0.5,
                               pred.reps = 1e2)

Fclst_raster.plot <- prevmaprasterplotter(Fclst_raster$pred,
                                          alpha = 1, smoothfct = 8)
load("data/map_bases/space_mips_maps_bases.rda")
drccites <- readr::read_csv("data/map_bases/DRC_city_coordinates.csv") %>%
  dplyr::filter(population > 350000)

jpeg("~/Desktop/temp.jpg", width = 11, height = 8, units = "in", res = 500)
Fclst_raster.plot +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1) +
  prettybasemap_nodrc_nonorth_dark +
  geom_point(data = drccites, aes(x = longnum, y=latnum), alpha = 0.5) +
  geom_text(data = drccites, aes(label = city, x = longnum, y=latnum),
            hjust = 0.5, vjust = 0.5, nudge_y = 0.25, fontface = "bold", alpha = 0.8)
graphics.off()





svglite::svglite("~/Desktop/temp.jpg", width = 8, height = 8)
Fclst_raster.plot +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1) +
  prettybasemap_nodrc_nonorth_dark +
  geom_point(data = drccites, aes(x = longnum, y=latnum), alpha = 0.5) +
  geom_text(data = drccites, aes(label = city, x = longnum, y=latnum),
            hjust = 0.5, vjust = 0.5, nudge_y = 0.25, fontface = "bold", alpha = 0.8)
graphics.off()

pfinc <- raster::raster("data/derived_data/MAPrasters/pfincidence.grd")
jpeg("~/Desktop/temp2.jpg", width = 11, height = 8, units = "in", res = 500)
ggplot() +
  ggspatial::layer_spatial(data = pfinc, aes(fill = stat(band1))) +
  scale_fill_distiller("Prevalence", type = "div", palette = "RdYlBu") +
  prettybasemap_nodrc_nonorth_dark +
  geom_point(data = drccites, aes(x = longnum, y=latnum), alpha = 0.5) +
  geom_text(data = drccites, aes(label = city, x = longnum, y=latnum),
            hjust = 0.5, vjust = 0.5, nudge_y = 0.25, fontface = "bold", alpha = 0.8)
graphics.off()

svglite::svglite("~/Desktop/temp2.svg", height = 8, width = 8)
ggplot() +
  ggspatial::layer_spatial(data = pfinc, aes(fill = stat(band1))) +
  scale_fill_distiller("Prevalence", type = "div", palette = "RdYlBu") +
  prettybasemap_nodrc_nonorth_dark +
  geom_point(data = drccites, aes(x = longnum, y=latnum), alpha = 0.5) +
  geom_text(data = drccites, aes(label = city, x = longnum, y=latnum),
            hjust = 0.5, vjust = 0.5, nudge_y = 0.25, fontface = "bold", alpha = 0.8)
graphics.off()

# sanity
