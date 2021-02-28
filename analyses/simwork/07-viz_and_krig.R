## .................................................................................
## Purpose: Visualize and Krig the estimates
##
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
library(fields)
library(ggspatial)
source("R/space_functions.R")

#............................................................
# read in results
#...........................................................
simmap <- readRDS("results/sim_inbreed_ests/min_cost_inbreedingresults/sim_min_cost_inbreedingresults.RDS")
locats <- readRDS("data/sim_data/lattice_model.rds") %>%
  dplyr::select(-c("migration"))

#............................................................
# run point plot functions
#...........................................................
simmap <- simmap %>%
  dplyr::mutate(pointplot = purrr::map(discentret, point_plotter, locats = locats),
                pointplot = purrr::map(pointplot, function(x){
                  x +
                    theme_void() +
                    scale_color_viridis_c("DISC", option = "cividis")
                }))

#............................................................
# run raster plot functions
#...........................................................
simmap <- simmap %>%
  dplyr::mutate(rasterpreds = purrr::map(discentret, get_krigged, locats = locats, res = 1e3),
                rasterplot = purrr::map(rasterpreds, raster_plotter),
                rasterplot = purrr::map(rasterplot, function(x){
                  x +
                    scale_color_viridis_c("Pred. DISC", option = "cividis") +
                    theme_void()
                }))



simdat <- readRDS("data/sim_data/sim_gengeodat.rds") %>%
  dplyr::select(c("q", "pairPlotObj"))

simmap <- simmap %>% dplyr::left_join(., simdat, by = "q")


simmap$pairPlotObj[[1]]
View(simmap$rasterpreds[[1]])
simmap$pointplot[[1]]
simmap$rasterplot[[1]]


simmap$pairPlotObj[[3]]
simmap$pointplot[[3]]
simmap$rasterplot[[3]]


simmap$pointplot[[2]]
simmap$pairPlotObj[[2]]
simmap$rasterplot[[2]]

simmap$pairPlotObj[[4]]
simmap$pointplot[[4]]
simmap$rasterplot[[4]]

