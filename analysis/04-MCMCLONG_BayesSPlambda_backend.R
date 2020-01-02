#..............................................................
# Purpose of this script is to fit the various models
# from CarBayes for our IBD Spatial Covariates.
# We use this "backend" script to interact with a server
#..............................................................
library(tidyverse)
library(CARBayes)
library(raster)
source("R/basics.R")
source("R/pairwise_helpers.R")
set.seed(48)

drcsmpls <- readRDS("data/distance_data/drcsmpls_foruse.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("id", "hv001")) %>%
  dplyr::rename(name = id)

mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(name = id) %>%
  dplyr::select(c("name", "country", "hv001", "adm1name", "longnum", "latnum")) %>%
  dplyr::filter(name %in% drcsmpls$name)

# drc prov for plot
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))


#....................................................................................
# Import Genetic Data
#....................................................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD_polarized_biallelic_processed.long.rds")

# log2 IBD measure
ibD <- ibD %>%
  dplyr::mutate(malecotf_gens = -log2(malecotf),
                malecotf_gens_inv = 1/malecotf_gens)
#..............................................................
# mk mtdt
#..............................................................
mtdt.clst <- mtdt %>%
  dplyr::select(c("name", "hv001")) %>%
  dplyr::rename(smpl1 = name)

ibD <- dplyr::left_join(ibD, mtdt.clst, by = "smpl1")

# rename
mtdt.clst <- mtdt.clst %>%
  dplyr::rename(smpl2 = smpl1)
# rejoin
ibD <- dplyr::left_join(ibD, mtdt.clst, by = "smpl2")


#..............................................................
# Import geographic distance and combine with genetic
#..............................................................
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")

ibDdist <- long_distance_matrix_join(x=ibD, y=distancematrix.cluster,
                                     by = c("hv001.x", "hv001.y")) %>%
  dplyr::mutate(gcdistance = gcdistance/1e3,
                roaddistance = roaddistance/1e3,
                riverdist = riverdist/1e3,
                riverdist = as.numeric(riverdist)) # remove units


# take care of cluster diagnonals
ibDdist$gcdistance[ibDdist$hv001.x == ibDdist$hv001.y] <- 0
ibDdist$roaddistance[ibDdist$hv001.x == ibDdist$hv001.y] <- 0
ibDdist$riverdist[ibDdist$hv001.x == ibDdist$hv001.y] <- 0

ibDdist.long <- ibDdist %>%
  dplyr::select(-c(dplyr::starts_with("hv001"))) %>%
  expand_distance_matrix(.)

p1 <- mtdt %>%
  dplyr::select(c("name", "adm1name")) %>%
  dplyr::rename(smpl1 = name)

p2 <- mtdt %>%
  dplyr::select(c("name", "adm1name")) %>%
  dplyr::rename(smpl2 = name)

ibDdist.long.mtdt <- ibDdist.long %>%
  dplyr::left_join(., y = p1, by = "smpl1") %>%
  dplyr::left_join(., y = p2, by = "smpl2")




#..............................................................
# Get within IBD for the Prov
#..............................................................
withinprovIBD <- ibDdist.long.mtdt %>%
  dplyr::filter(adm1name.x == adm1name.y) %>%
  dplyr::group_by(adm1name.x) %>%
  dplyr::summarise(
    meangens = mean(malecotf)
  ) %>%
  dplyr::rename(adm1name = adm1name.x)




#..............................................................
# Import Covariates for Province
#..............................................................
provCovar.raw <- readRDS("data/derived_data/covar_rasterstack_provlocations_raw.RDS")

# transform covars
provCovar <- provCovar.raw %>%
  dplyr::mutate(
    prev = my.scale(logit(prev, tol = 0.1)),
    precip = my.scale(precip),
    temp = my.scale(temp),
    elev = my.scale(elev),
    crops = my.scale(logit(crops, tol = 0.1)),
#    actuse = my.scale(logit(actuse, tol = 0.1)), actuse doesn't have any real variation in DRC
    netuse = my.scale(logit(netuse, tol = 0.1)),
    housing = my.scale(logit(housing, tol = 0.1))
  )


#..............................................................
# Combine outcome and covariate data
#..............................................................
# housekeeping
riskvar <- colnames(provCovar)[!colnames(provCovar) %in% "adm1name"]
# model framework
mod.IBD.provCovar <- dplyr::left_join(withinprovIBD, provCovar,
                                       by = "adm1name")




#-------------------------------------------------------------------------
# Conditional Autoregressive Spatial Model
#-------------------------------------------------------------------------
#..............................................................
# Make Adjacency Matrix
# by space
#..............................................................
#W.nb <- spdep::poly2nb(sf::as_Spatial(DRCprov), row.names = DRCprov$adm1name)
#W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)

make_symm_mat <- function(x){
  # now make it a symm matrix
  ret <- rbind( matrix(NA, nrow = 1, ncol = ncol(x)), as.matrix(x) )
  ret <- cbind(ret, matrix(NA, nrow = nrow(ret), ncol = 1) )
  diag(ret) <- 0
  ret[upper.tri(ret)] <-  t(ret)[upper.tri(ret)]
  return(ret)
}
#..............................................................
# get distance matrices
#..............................................................
# gc
prov.gcdist <- readRDS("data/distance_data/greater_circle_distance_forprovinces.rds")
W.gcdist <- prov.gcdist %>%
  tidyr::spread(., key = "item2", value = "gcdistance") %>%
  dplyr::select(-c("item1"))
# now make it a symm matrix
W.gcdist <- make_symm_mat(W.gcdist)/1e3 # meters to kilometers, to help with variance

# road
prov.roaddist <- readRDS("data/distance_data/prov_road_distmeters_long.rds")
W.roaddist <- prov.roaddist %>%
  tidyr::spread(., key = "item2", value = "roaddistance") %>%
  dplyr::select(-c("item1"))
# now make it a symm matrix
W.roaddist <- make_symm_mat(W.roaddist)/1e3 # meters to kilometers, to help with variance

# river
prov.riverdist <- readRDS("data/distance_data/river_distance_forprovinces.rds")
W.riverdist <- prov.riverdist %>%
  dplyr::select(c("dhsprovfrom", "dhsprovto", "riverdist")) %>%
  dplyr::rename(item1 = dhsprovfrom,
                item2 = dhsprovto) %>%
  tidyr::spread(., key = "item2", value = "riverdist") %>%
  dplyr::select(-c("item1"))
# now make it a symm matrix
W.riverdist <- make_symm_mat(W.riverdist)/1e3 # meters to kilometers, to help with variance


#..............................................................
# Run Leroux model in parallel on slurm
#..............................................................
wrap_S.CARleroux <- function(distcat, outcome, formula, W, data, burnin, n.sample){

  #..............................................................
  # fit model
  #..............................................................
  formvec <- paste(deparse(formula), collapse = "")
  betacount <- stringr::str_count(formvec, "\\+") + 1 # plus one for intercept
  prior.var.betavec <- rep(5e4, betacount) # note prior setting here

  ret <- CARBayes::S.CARleroux(formula = formula,
                               family = "gaussian",
                               W = W,
                               data = data,
                               burnin = burnin,
                               prior.var.beta = prior.var.betavec,
                               prior.tau2 = c(1, 0.01),
                               n.sample = n.sample)
  # note, rho by default is NULL which allows it to be estimated in the model


  #-------------------------------------------------------------------------
  # MCMC Diagnostics
  #-------------------------------------------------------------------------
  ret <- tibble::tibble(MCMC = list(ret))
  ret$mcmc.modsum <- purrr::map(ret$MCMC, print)  # note, print is overloaded here
  ret$summresults <- purrr::map(ret$mcmc.modsum, "summary.results")
  ret$summresults <- purrr::map(ret$summresults, function(x){
    pars <- rownames(x)
    ret <- cbind.data.frame(pars = pars, as.data.frame(x))
    return(ret)
  })

  #..............................................................
  # DIC for fits
  #..............................................................
  ret$modfit <- purrr::map(ret$mcmc.modsum, "modelfit")
  ret$DIC <- purrr::map(ret$modfit, "DIC")

  #..............................................................
  # Keep smaller piece
  #..............................................................
  ret <- ret %>%
    dplyr::select(c("summresults", "DIC")) %>%
    tidyr::unnest(cols = c("summresults", "DIC"))

  return(ret)

}

#..............................................................
# Spatial random effects Models
#..............................................................
mod.framework.sp <- tibble(distcat = c("gc", "road", "river"),
                           formula = list(as.formula("meangens ~ netuse + 1"),
                                          as.formula("meangens ~ netuse + 1"),
                                          as.formula("meangens ~ 1")),
                           burnin = 1e4,
                           n.sample = 1e7 + 1e4,
                           W = list(W.gcdist, W.roaddist, W.riverdist),
                           data = list(mod.IBD.provCovar))

mod.framework.sp$results <- purrr::pmap(mod.framework.sp, wrap_S.CARleroux)
dir.create("results/carbayes_longchains")
saveRDS(mod.framework.sp, "results/carbayes_longchains/longchains_models.RDS")

