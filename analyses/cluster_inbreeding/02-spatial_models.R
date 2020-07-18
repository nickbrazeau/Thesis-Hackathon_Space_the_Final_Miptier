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

#............................................................
# function for spatial model
#...........................................................
make_spat_mod <- function(clst_inbdset, DRCprov, clsts, covar = "1", kappa = 0.5) {
  #......................
  # process
  #......................
  clst_inbdset <- clst_inbdset %>%
    dplyr::filter(param != "m") %>%
    dplyr::rename(hv001 = param,
                  Finbd = est)
  clst_inbdset <- dplyr::left_join(clst_inbdset, clsts, by = "hv001") %>%
    dplyr::mutate(Finbd_logit = logit(as.numeric(Finbd)))
  #......................
  # run internal prevmap functions
  #......................
  poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14))
  grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
  colnames(grid.pred) <- c("long","lat")

  Fclst_raster <- fit_pred_spMLE(data = clst_inbdset,
                                 outcome = "Finbd_logit", covar = covar,
                                 long_var = "longnum", lat_var = "latnum",
                                 grid.pred = grid.pred, kappa = kappa,
                                 pred.reps = 1e3)

  Fclst_raster_plot <- prevmaprasterplotter(Fclst_raster$pred,
                                            alpha = 1, smoothfct = 8)

  Fclst_raster_plot_obj <- Fclst_raster_plot +
    scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1) +
    prettybasemap_nodrc_nonorth_dark +
    geom_point(data = drccites, aes(x = longnum, y=latnum), alpha = 0.5) +
    geom_text(data = drccites, aes(label = city, x = longnum, y=latnum),
              hjust = 0.5, vjust = 0.5, nudge_y = 0.25, fontface = "bold",
              size = 3,
              alpha = 0.8)

  return(Fclst_raster_plot_obj)
}

#..............................................................
# read in data
#..............................................................
clst_inbd <- readRDS("results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = inbreed_ests)
clst_inbd.list <- split(clst_inbd, factor(clst_inbd$spacetype))

#......................
# geo
#......................
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))
clsts <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("hv001", "latnum", "longnum")) %>%
  dplyr::mutate(hv001 = as.character(hv001)) %>%
  dplyr::filter(!duplicated(.))
load("data/map_bases/space_mips_maps_bases.rda")
drccites <- readr::read_csv("data/map_bases/DRC_city_coordinates.csv") %>%
  dplyr::filter(population > 350000)



#..............................................................
# run functions
#..............................................................
clst_inbd.plots <- lapply(clst_inbd.list, make_spat_mod, clsts = clsts, DRCprov = DRCprov)

jpeg("~/Desktop/inbred_fig.jpg", width = 11, height = 8, units = 'in', res = 600)
cowplot::plot_grid(clst_inbd.plots[[1]], clst_inbd.plots[[2]],
                   clst_inbd.plots[[3]], clst_inbd.plots[[4]],
                   nrow = 2, labels = c("(A)", "(B)", "(C)", "(D)"))
graphics.off()








# sanity
