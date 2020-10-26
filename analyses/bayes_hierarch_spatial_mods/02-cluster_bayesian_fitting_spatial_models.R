#......................................................................................................
## Purpose: Base Spatial Modeling of Cluster Inbreeding Coeffs
##
## Author: Nicholas F. Brazeau
##
## Date: April 08 2020
#......................................................................................................
source("R/basics.R")
library(tidyverse)

#......................................................................................................
####  Point Process-Based Modeling ####
#......................................................................................................
library(PrevMap)
library(rgeos)
library(gstat)
#......................
# bring in data
#......................
# covariate data
clst_covars <- readRDS("data/derived_data/covar_rasterstack_samplinglocations_scaled.RDS") %>%
  dplyr::mutate(hv001 = as.character(hv001))

# load pretty map aesthetics
load("data/map_bases/space_mips_maps_bases.rda")

# spatial data
ge <- readRDS("data/derived_data/spacemips_GE.rds") %>%
  dplyr::mutate(hv001 = as.character(hv001))
ge <- sf::st_as_sf(ge)
sf::st_crs(ge) <-  "+proj=longlat +datum=WGS84 +no_defs"


# inbreeding data
clst_inbd <- readRDS("results/clust_inbd_results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS") %>%
  dplyr::filter(spacetype != "migrate") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = inbreed_ests) %>%
  dplyr::mutate(logitest = logit(est)) %>%
  dplyr::ungroup(.) %>%
  dplyr::rename(hv001 = param)

# bring together
clst_spatdat <- dplyr::left_join(clst_inbd, clst_covars, by = "hv001") %>%
  dplyr::left_join(., ge, by = "hv001") %>%
  dplyr::filter(hv001 != "m") %>%
  dplyr::group_by(spacetype) %>%
  tidyr::nest(.)

#............................................................
# Modelling
#...........................................................
#......................
# going to assume that kappa is 0.5 and follows and exponential model -- based
# on the isolation by distance framework
# find appropriate start parameters for tau, sigma, phi
#......................
find_spat_covars <- clst_spatdat %>%
  dplyr::mutate(geoRobj = purrr::map(data, geoR::as.geodata, # make GeoR Object
                                     coords.col = c("longnum", "latnum"),
                                     data.col = "logitest",
                                     covar.col = c("incidence", "urban")),
                # exponential variogram and fit
                variog_exp = purrr::map(geoRobj, geoR::variog,
                                        uvec =  c(seq(0, 2, by = 0.1))),
                variofit_exp = purrr::map(variog_exp, geoR::variofit,
                                          cov.model = "exp",
                                          ini = c(2, 0.1),
                                          nugget = 0,
                                          fix.nugget = F,
                                          fix.kappa = T,
                                          kappa = 0.5),
                # matern variogram and fit
                variog_matern = purrr::map(geoRobj, geoR::variog,
                                           uvec =  c(seq(0, 2, by = 0.1))),
                variofit_matern = purrr::map(variog_matern, geoR::variofit,
                                             cov.model = "matern",
                                             ini = c(2, 0.1),
                                             nugget = 0,
                                             fix.nugget = F,
                                             fix.kappa = T,
                                             kappa  = 0.5))


find_spat_covars$variofit_exp[[1]]
find_spat_covars$variofit_exp[[2]]
find_spat_covars$variofit_matern[[1]]
find_spat_covars$variofit_matern[[2]]


#......................
# Model framework
#......................
clst_mod_framework <- tibble(spacetype = c(rep("gcdist", 2), rep("roaddist", 2)),
                             name = rep(c("intercept",
                                          "interaction"), 2),
                             formula = rep(c("logitest ~ 1",
                                             "logitest ~ incidence + urban + incidence*urban"), 2))


# bring in data
clst_mod_framework <- clst_mod_framework %>%
  dplyr::left_join(., clst_spatdat) %>%
  dplyr::mutate(fitlm = purrr::map2(formula, data, function(x, y){lm(as.formula(x),
                                                                     data = y)}))
# make prior information
clst_mod_framework <- clst_mod_framework %>%
  dplyr::mutate(mypriors = purrr::map(formula, function(x){
    formvec <- paste(x, collapse = "")
    betacount <- stringr::str_count(formvec, "\\+") + 2 # need intercept and betas
    betacount <- ifelse(grepl("1", formvec), 1, betacount) # corner case of just intercept
    if (betacount > 1) {
      betamean <- rep(0, betacount)
      covarsmat <- matrix(0, ncol = betacount, nrow=betacount)
      diag(covarsmat) <- 1 # identity matrix
    } else {
      betamean <- 0
      covarsmat <- 1
    }
    mypriors <- PrevMap::control.prior(beta.mean = betamean,
                                       beta.covar = covarsmat,
                                       uniform.nugget = c(0,1), # this is tau2
                                       uniform.phi = c(6,18),
                                       log.normal.sigma = c(-1.3, 1))
    return(mypriors)
  }))

# make directions for mcmc
clst_mod_framework <- clst_mod_framework %>%
  dplyr::mutate(mydirections = purrr::map2(formula, fitlm, function(x, y){
    formvec <- paste(x, collapse = "")
    betacount <- stringr::str_count(formvec, "\\+") + 2 # need intercept and betas
    betacount <- ifelse(grepl("1", formvec), 1, betacount) # corner case of just intercept

    mydirections <- PrevMap::control.mcmc.Bayes(burnin = 1e3,
                                                n.sim = 1e3+1e3,
                                                thin = 10,
                                                L.S.lim = c(5,50),
                                                epsilon.S.lim = c(0.01, 0.1),
                                                start.nugget = 0.05,
                                                start.sigma2 = 1,
                                                start.beta = rep(0, betacount),
                                                start.phi = 13,
                                                start.S = predict(y),
                                                linear.model = TRUE)
    return(mydirections)
  }))

#......................
# Make a knot bounding box
#......................
drcprov <- readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")
bb <- as(raster::extent(11, 32,-15, 7), "SpatialPolygons")
sp::proj4string(bb) <- "+proj=longlat +datum=WGS84 +no_defs"
ggplot() +
  geom_sf(data = drcprov) +
  geom_sf(data = sf::st_as_sf(bb), color = "red", alpha = 0.2)
# want bounding box outside borders for boundary effects
knotsbb <- expand.grid(seq(11, -15, length.out = 10), seq(32, 7, length.out = 10))

#......................
# Make a wrapper for PrevMap fitting
#......................
fit_bayesmap_wrapper <- function(formula,
                                 data,
                                 mypriors,
                                 mydirections,
                                 kappa,
                                 knots){

  ret <- PrevMap::linear.model.Bayes(
    formula = as.formula(formula),
    coords = as.formula("~ longnum + latnum"),
    data = data,
    control.prior = mypriors,
    control.mcmc = mydirections,
    kappa = kappa,
    low.rank = TRUE,
    knots = knots)
  return(ret)
}

# run
clst_mod_framework$PrevMapFit <- purrr::pmap(clst_mod_framework[, c("formula",
                                                                    "data",
                                                                    "mypriors",
                                                                    "mydirections")],
                                             fit_bayesmap_wrapper,
                                             kappa = 0.5, # exponential model
                                             knots = knotsbb)

# save out
saveRDS(clst_mod_framework, "results/clust_inbd_results/final_clstb_maps/clst_mod_framework_prevmapfits.RDS")


#............................................................
# Predictions
#...........................................................
clst_mod_preds <- clst_mod_framework %>%
  dplyr::select(c("spacetype", "name", "PrevMapFit"))

#......................
# do this hack to tell the model the formula it call
#......................
clst_mod_preds <- clst_mod_preds %>%
  dplyr::mutate(PrevMapFit = ifelse(name == "intercept",
                                    purrr::map(PrevMapFit, function(x){
                                      x$call$formula = as.formula("logitest ~ 1")
                                      return(x)
                                    }),
                                    purrr::map(PrevMapFit, function(x){
                                      x$call$formula = as.formula("logitest ~ incidence + urban + incidence*urban")
                                      return(x)
                                    })))

#...............................
# sample coordinates for prediction surface
#...............................
# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14))
grid.pred.coords <- splancs::gridpts(poly, xs=0.05, ys=0.05)
colnames(grid.pred.coords) <- c("longnum","latnum")
grid.pred.coords.df <- as.data.frame(grid.pred.coords)
# Raster surfaces for risk factors
covars <- c("incidence",
            "urban", "incidence:urban")
covar.rasterstack.derived <- readRDS("data/derived_data/covar_rasterstack_derived.RDS")
inc_urban <- raster::overlay(covar.rasterstack.derived, fun = function(x,y){x*y})
names(inc_urban) <- "incidence:urban"
predcovars <- raster::stack(covar.rasterstack.derived, inc_urban)
# pred df
pred.df <- raster::extract(
  x = predcovars,
  y = sf::as_Spatial(
    sf::st_as_sf(grid.pred.coords.df, coords = c("longnum", "latnum"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs")),
  buffer = 6000,
  fun = mean,
  na.rm = T,
  sp = F)

pred.df <- as.data.frame(pred.df)
pred.df <- cbind.data.frame(grid.pred.coords.df, pred.df)

# down sample
covar.rstr.pred <- raster::rasterFromXYZ(pred.df,
                                         res = c(0.05, 0.05),
                                         crs = "+proj=longlat +datum=WGS84 +no_defs")

covar.rstr.pred.downsmpl <- raster::sampleRandom(covar.rstr.pred,
                                                 size = 2e4,
                                                 na.rm = T,
                                                 xy = T,
                                                 sp = F)

colnames(covar.rstr.pred.downsmpl) <- c("longnum", "latnum",
                                        "incidence", "urban", "incidence:urban")
colnames(covar.rstr.pred.downsmpl)[5] <- "incidence:urban"

#...............................
# Setup predictions
#...............................
# set up grid.pred
clst_mod_preds$grid.pred <- list(covar.rstr.pred.downsmpl[,c("longnum", "latnum")])

# set up predictors
clst_mod_preds <- clst_mod_preds %>%
  dplyr::mutate(predictors = ifelse(name == "intercept", list(NULL),
                                    list(as.data.frame(
                                      covar.rstr.pred.downsmpl[,c("incidence",
                                                                  "urban",
                                                                  "incidence:urban")])))
  )


# make wrapper
pred_PrevMap_bayes_wrapper <- function(PrevMapFit, grid.pred, predictors){
  ret <- PrevMap::spatial.pred.linear.Bayes(object = PrevMapFit,
                                            grid.pred = grid.pred,
                                            predictors = predictors,
                                            type = "marginal",
                                            scale.predictions = "prevalence",
                                            quantiles = NULL,
                                            standard.error = T,
                                            thresholds = NULL)
  return(ret)

}

# get predictions
clst_mod_preds$PrevMapPredictions <- purrr::pmap(clst_mod_preds[,c("PrevMapFit", "grid.pred", "predictors")],
                                                 pred_PrevMap_bayes_wrapper)


#......................
# tidy up predictions
#......................
clst_mod_preds <- clst_mod_preds %>%
  dplyr::mutate(PrevMap_postmeans = purrr::map(PrevMapPredictions, function(x){
    clstmodel.fitted <- expit(x$samples)
    clstmodel.coords <- x$grid.pred
    clstmodel.fitted.postmeans <- apply(clstmodel.fitted, 2, mean)
    clstmodel.coords.means <- cbind.data.frame(clstmodel.coords,
                                               fitted.postmean = clstmodel.fitted.postmeans)
    return(clstmodel.coords.means)
  }))

# save out
saveRDS(clst_mod_preds, "results/clust_inbd_results/final_clstb_maps/clst_mod_framework_prevmappreds.RDS")


#......................
# do local interpolation to fill in rest of points
#......................
# for masking
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
# for interpolating
drcrstr <- raster::rasterFromXYZ(cbind(grid.pred.coords[,1],
                                       grid.pred.coords[,2],
                                       NA),
                                 crs="+proj=longlat +datum=WGS84 +no_defs")

do_local_interpolate_toPlotObj <- function(PrevMap_postmeans) {
  # local interpolate
  colnames(PrevMap_postmeans)[1:2] <- c("x", "y")
  PrevMap_postmeans.idwod <- gstat::gstat(id = "postmeans", formula = fitted.postmean ~ 1,
                                          locations = ~x + y,
                                          data = PrevMap_postmeans,
                                          set=list(idp = 2))
  PrevMap_postmeans.idwod <- raster::interpolate(drcrstr, PrevMap_postmeans.idwod)
  PrevMap_postmeans.idwod <- raster::mask(PrevMap_postmeans.idwod, DRC)

  # plot out
  plotObj <- ggplot() +
    ggspatial::layer_spatial(data = PrevMap_postmeans.idwod,
                             aes(fill = stat(band1))) +
    scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1, na.value = NA) +
    prettybasemap_nodrc_nonorth_dark +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0),
      legend.text = element_text(face = "bold", size = 11),
      plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))
  return(plotObj)
}

# do work
clst_mod_preds$plotObj <- purrr::map(clst_mod_preds$PrevMap_postmeans,
                                     do_local_interpolate_toPlotObj)

# save
saveRDS(clst_mod_preds$plotObj[[1]],
        file = "results/clust_inbd_results/final_clstb_maps/intercept_gc_clst_mod_framework_prevmapPlot.RDS")
saveRDS(clst_mod_preds$plotObj[[2]],
        file = "results/clust_inbd_results/final_clstb_maps/interaction_gc_clst_mod_framework_prevmapPlot.RDS")
saveRDS(clst_mod_preds$plotObj[[3]],
        file = "results/clust_inbd_results/final_clstb_maps/intercept_road_clst_mod_framework_prevmapPlot.RDS")
saveRDS(clst_mod_preds$plotObj[[4]],
        file = "results/clust_inbd_results/final_clstb_maps/interaction_road_clst_mod_framework_prevmapPlot.RDS")




