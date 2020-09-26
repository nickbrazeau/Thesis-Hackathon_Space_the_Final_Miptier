#########################################################################
# Purpose: Base Spatial Modeling of Cluster Inbreeding Coeffs
#
# Author: Nicholas F. Brazeau
#
# Date: April 08 2020
#########################################################################
source("R/basics.R")
library(tidyverse)
library(gstat)

# load pretty map aesthetics
load("data/map_bases/space_mips_maps_bases.rda")
# add cities for context
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))
clsts <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("hv001", "latnum", "longnum")) %>%
  dplyr::mutate(hv001 = as.character(hv001)) %>%
  dplyr::filter(!duplicated(.))

#............................................................
# functions for spatial model
#...........................................................
#' @title Raw Map Point Process
make_spat_raw_map <- function(clst_inbdset, DRCprov, clsts, covar = "1", kappa = 0.5) {
  #......................
  # process
  #......................
  clst_inbdset <- clst_inbdset %>%
    dplyr::filter(param != "m") %>%
    dplyr::rename(hv001 = param,
                  Finbd = est)
  clst_inbdset <- dplyr::left_join(clst_inbdset, clsts, by = "hv001")

  Fclst_point_plot_obj <- clst_inbdset %>%
    ggplot() +
    geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
    geom_point(aes(x = longnum, y = latnum, color = Finbd)) +
    scale_color_viridis_c("Inbreeding", option="plasma", direction = 1)

  return(Fclst_point_plot_obj)
}


#' @title Inverse Distance Weighting Kriging
interpolate_spat_idw_mod <- function(clst_inbdset, DRCprov, clsts, idw = 2) {
  #......................
  # process
  #......................
  clst_inbdset <- clst_inbdset %>%
    dplyr::filter(param != "m") %>%
    dplyr::rename(hv001 = param,
                  Finbd = est)
  clst_inbdset <- dplyr::left_join(clst_inbdset, clsts, by = "hv001")
  #......................
  # get internal pieces needed for interpolation
  #......................
  poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14))
  grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
  colnames(grid.pred) <- c("longnum","latnum")

  # need this for bounding box
  DRC <- sf::as_Spatial(osmdata::getbb("Democratic Republic of the Congo",
                                       featuretype = "country",
                                       format_out = 'sf_polygon'))
  drcrstr <- raster::rasterFromXYZ(cbind(grid.pred[,1],
                                         grid.pred[,2],
                                         NA),
                                   crs="+proj=longlat +datum=WGS84 +no_defs")

  #..............................................................
  # interpolate
  #..............................................................

  idwmod <- gstat::gstat(formula = Finbd ~ 1,
                         locations = ~longnum + latnum,
                         data = clst_inbdset,
                         set=list(idp = idw))

  ret <- raster::interpolate(drcrstr, idwmod, xyNames = c("longnum", "latnum"))
  ret <- raster::mask(ret, DRC)

  # plot object
  Fclst_raster_plot_obj <-  ggplot() +
    ggspatial::layer_spatial(data = ret,
                             aes(fill = stat(band1)),
                             alpha = alpha) +
    scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1) +
    prettybasemap_nodrc_nonorth_dark +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0),
      legend.text = element_text(face = "bold", size = 11),
      plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))
  return(Fclst_raster_plot_obj)
}

#..............................................................
# read in data cluster
#..............................................................
# cluster names
clsts <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("hv001", "longnum", "latnum")) %>%
  dplyr::filter(!duplicated(.)) %>%
  tibble::as_tibble(.) %>%
  dplyr::mutate(hv001 = as.character(hv001))

# inbreeding data
clst_inbd <- readRDS("results/clust_inbd_results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS") %>%
  dplyr::filter(spacetype != "migrate") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = inbreed_ests)
clst_inbd.list <- split(clst_inbd, factor(clst_inbd$spacetype))



#..............................................................
# run functions for clusters
#..............................................................
#......................
# raw
#......................
clst_inbd_point.plots <- lapply(clst_inbd.list, make_spat_raw_map, clsts = clsts, DRCprov = DRCprov)

#......................
# inverses distance weighting
#......................
clst_inbd_rstr.plots <- lapply(clst_inbd.list, interpolate_spat_idw_mod, clsts = clsts, DRCprov = DRCprov)

# out
dir.create("results/clust_inbd_results/final_maps/", recursive = T)
saveRDS(clst_inbd_point.plots, file = "results/clust_inbd_results/final_clstb_maps/point_clst_inbd_plots.RDS")
saveRDS(clst_inbd_rstr.plots, file = "results/clust_inbd_results/final_clstb_maps/raster_clst_inbd_plots.RDS")




#..............................................................
# read in data for voroni tesselated territories
#..............................................................
vrdf <- readRDS("data/distance_data/voroni_base.RDS") %>%
  dplyr::rename(param = IPUMSID) %>%
  dplyr::mutate(param = as.character(param)) %>%
  rename(geometry = x)
sf::st_crs(vrdf) <-  "+proj=longlat +datum=WGS84 +no_defs"

# inbreeding data
prov_inbd <- readRDS("results/clust_inbd_results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS") %>%
  dplyr::filter(spacetype == "migrate") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = inbreed_ests) %>%
  dplyr::ungroup(.)
#......................
# plot out
#......................
prov_inbd_vrdf <- dplyr::left_join(vrdf, prov_inbd)
prov_inbd_vrdf_plotObj <- ggplot() +
  geom_sf(data = prov_inbd_vrdf, aes(fill = est)) +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1) +
  prettybasemap_nodrc_nonorth_dark +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0),
    legend.text = element_text(face = "bold", size = 11),
    plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))

# save
saveRDS(prov_inbd_vrdf, file = "results/clust_inbd_results/final_maps/prov_inbd_vrdf_terrmap.RDS")





# sanity
