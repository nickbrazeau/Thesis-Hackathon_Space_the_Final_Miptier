#......................................................................................................
## Purpose: Base Spatial Modeling of Cluster Inbreeding Coeffs
##
## Author: Nicholas F. Brazeau
##
## Date: April 08 2020
#......................................................................................................
source("R/basics.R")
library(tidyverse)
load("data/map_bases/space_mips_maps_bases.rda")

#............................................................
# raw maps
#...........................................................
provraw <- readRDS("results/clust_inbd_results/final_prov_maps/bayesian_posterior_prov_inbd_intercept_plot.RDS") +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1, na.value = NA) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(face = "bold", size = 11, angle = 45, hjust = 0.5, vjust = 0),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))

rawgcdist <- readRDS("results/clust_inbd_results/final_clstb_maps/intercept_gc_clst_mod_framework_prevmapPlot.RDS") +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1, na.value = NA) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(face = "bold", size = 11, angle = 45, hjust = 0.5, vjust = 0),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))

rawroaddist <- readRDS("results/clust_inbd_results/final_clstb_maps/intercept_road_clst_mod_framework_prevmapPlot.RDS") +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1, na.value = NA) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(face = "bold", size = 11, angle = 45, hjust = 0.5, vjust = 0),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))
#............................................................
# fitted maps
#...........................................................
provmodelplotObj <- readRDS("results/clust_inbd_results/final_prov_maps/bayesian_posterior_prov_inbd_interaction_plot.RDS")
provmodelplotObj <- provmodelplotObj + prettybasemap_nodrc_nonorth_dark +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1, na.value = NA) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(face = "bold", size = 11, angle = 45, hjust = 0.5, vjust = 0),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))

fitted_gc_prevmeans <- readRDS("results/clust_inbd_results/final_clstb_maps/interaction_gc_clst_mod_framework_prevmapPlot.RDS") +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1, na.value = NA) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(face = "bold", size = 11, angle = 45, hjust = 0.5, vjust = 0),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))

fitted_road_prevmeans <- readRDS("results/clust_inbd_results/final_clstb_maps/interaction_road_clst_mod_framework_prevmapPlot.RDS") +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1, na.value = NA) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(face = "bold", size = 11, angle = 45, hjust = 0.5, vjust = 0),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))

#............................................................
# bringtogether
#...........................................................
jpeg("results/figures/bayes_spatial_models_raw_fitted.jpg", width = 11, height = 8, units = "in", res = 600)
cowplot::plot_grid(rawgcdist, rawroaddist, provraw,
                   fitted_gc_prevmeans, fitted_road_prevmeans, provmodelplotObj,
                   nrow = 2, ncol = 3)
graphics.off()
