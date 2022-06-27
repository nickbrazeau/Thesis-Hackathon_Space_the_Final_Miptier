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
simmap <- readRDS("results/simulated_min_cost_inbreedingresults/sim_min_cost_inbreedingresults.RDS")
locats <- readRDS("data/sim_data/lattice_model.rds") %>%
  dplyr::select(-c("migration"))

#......................
# bring in realized ibd
#......................
simdat <- readRDS("data/sim_data/sim_gengeodat.rds")
simmap <- simmap %>% dplyr::left_join(., simdat, by = "q") %>%
  dplyr::mutate(lvl = factor(q, levels = c("mq_good", "mq_bad", "coi", "ne",
                                           "mtn", "oppcorner", "rift"))) %>%
  dplyr::arrange(lvl)



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
# do together discent and point plot
simmap$plotObj_realIBD_discent <- purrr::pmap(simmap[,c("gengeodat", "discentret")],
                                              plot_realIBD_discent, locats = locats)

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



#......................
# viz out
#......................
simmap$pairPlotObj[[1]]
simmap$pointplot[[1]]
simmap$plotObj_realIBD_discent[[1]]
simmap$rasterplot[[1]]

simmap$pairPlotObj[[2]]
simmap$pointplot[[2]]
simmap$plotObj_realIBD_discent[[2]]
simmap$rasterplot[[2]]

simmap$pairPlotObj[[3]]
simmap$pointplot[[3]]
simmap$plotObj_realIBD_discent[[3]]
simmap$rasterplot[[3]]

simmap$pairPlotObj[[4]]
simmap$pointplot[[4]]
simmap$plotObj_realIBD_discent[[4]]
simmap$rasterplot[[4]]

simmap$pairPlotObj[[5]]
simmap$pointplot[[5]]
simmap$plotObj_realIBD_discent[[5]]
simmap$rasterplot[[5]]

simmap$pairPlotObj[[6]]
simmap$pointplot[[6]]
simmap$plotObj_realIBD_discent[[6]]
simmap$rasterplot[[6]]

simmap$pairPlotObj[[7]]
simmap$pointplot[[7]]
simmap$plotObj_realIBD_discent[[7]]
simmap$rasterplot[[7]]


