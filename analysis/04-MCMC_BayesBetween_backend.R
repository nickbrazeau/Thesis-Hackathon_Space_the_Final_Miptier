#..............................................................
# Purpose of this script is to fit the various models
# from CarBayes for our IBD Spatial Covariates.
# We use this "backend" script to interact with a server
# OUTCOME: BETWEEN PROV
# Model: Y ~ B_{btwn} + B^T_{preds} + B_{\deltadist}
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

p1 <- mtdt %>%
  dplyr::select(c("name", "adm1name")) %>%
  dplyr::rename(smpl1 = name)

p2 <- mtdt %>%
  dplyr::select(c("name", "adm1name")) %>%
  dplyr::rename(smpl2 = name)



################################################################################
##############             Diagnostic Models         ###########################
################################################################################
#..............................................................
# Get between prov IBD
#..............................................................
# here we we don't need long format since we can parse on
# pairwise comparisons
ibDdist.mtdt <- ibDdist %>%
  dplyr::left_join(., y = p1, by = "smpl1") %>%
  dplyr::left_join(., y = p2, by = "smpl2")

provs <- sort( unique(c(as.character(ibDdist.mtdt$adm1name.x),
                  as.character(ibDdist.mtdt$adm1name.y))) )
# get between provs
provs.btwn <- t(combn(x = c(provs), m = 2))
# bind self
ibDdist.prov.between <- rbind.data.frame(provs.btwn,
                               unname(cbind(provs, provs)), stringsAsFactors = F) %>%
  dplyr::rename(prov.x = V1,
                prov.y = V2) %>%
  dplyr::mutate(btwn = paste0(prov.x,  "/", prov.y)) %>%
  dplyr::select(c("btwn", dplyr::starts_with("prov")))

#..............................................................
# Make Model Dataframe for Mapping Regression Models
#..............................................................
ibDdist.mtdt <- ibDdist.mtdt %>%
  dplyr::rename(riverdistance = riverdist) %>%
  dplyr::select(c("smpl1", "smpl2", "adm1name.x", "adm1name.y", "malecotf", "malecotf_gens", dplyr::ends_with("distance")))

ibDdist.prov.between$data <- purrr::pmap(ibDdist.prov.between[,c("prov.x", "prov.y")],
                                         function(prov.x, prov.y){
                                           if(prov.x == prov.y){ # same prov
                                             ret <- ibDdist.mtdt %>%
                                               dplyr::filter(adm1name.x == adm1name.y)
                                           } else { # for all others
                                             ret <- ibDdist.mtdt %>%
                                               dplyr::filter(
                                                 c(adm1name.x == prov.x & adm1name.y == prov.y |
                                                     adm1name.y == prov.x & adm1name.x == prov.y
                                                 )
                                               )
                                           }
                                           return(ret)
                                         }
                                         )

# make tibble for easier viz
ibDdist.prov.between <- tibble::as_tibble(ibDdist.prov.between)

#..............................................................
# Tidy and group
#..............................................................
ibDdist.prov.between$data <- purrr::map(ibDdist.prov.between$data,
                                        function(x){
                                          ret <- x %>%
                                            tidyr::gather(., key = "distcat", value = "geodist", 7:9)
                                          return(ret)
                                        })
# clean up df
ibDdist.prov.between <- ibDdist.prov.between %>%
  dplyr::select(-c("prov.x", "prov.y")) %>%
  tidyr::unnest(cols = "data") %>%
  dplyr::group_by(btwn, distcat) %>%
  tidyr::nest()

#..............................................................
# Get Mean IBD between Prov
#..............................................................
ibDdist.prov.between$F_between_mean <- purrr::map(ibDdist.prov.between$data, function(x){
  return(mean(x$malecotf))
})
ibDdist.prov.between <- ibDdist.prov.between%>%
  dplyr::ungroup(.)

ibDdist.prov.between <- ibDdist.prov.between[,c("btwn", "distcat", "F_between_mean")] %>%
  tidyr::unnest(cols = "F_between_mean")

#..............................................................
# Import Covariates for Province
#..............................................................
provCovar.raw <- readRDS("data/derived_data/covar_rasterstack_provlocations_raw.RDS")
# here we want differences by prov
ibDdist.prov.between$predMat <- purrr::map(ibDdist.prov.between$btwn, function(x){
  prov1 <- stringr::str_split_fixed(x, "/", n=2)[,1]
  prov2 <- stringr::str_split_fixed(x, "/", n=2)[,2]

  prov1 <- provCovar.raw %>%
    dplyr::filter(adm1name == prov1) %>%
    dplyr::select(-c("adm1name"))

  prov2 <- provCovar.raw %>%
    dplyr::filter(adm1name == prov2) %>%
    dplyr::select(-c("adm1name"))

  ret <- (prov1 - prov2)^2
  return(ret)

})


#..............................................................
# Get distance covariate
#..............................................................
# distances
prov.gcdist <- readRDS("data/distance_data/greater_circle_distance_forprovinces.rds")
prov.roaddist <- readRDS("data/distance_data/prov_road_distmeters_long.rds")
prov.riverdist <- readRDS("data/distance_data/river_distance_forprovinces.rds") %>%
  dplyr::select(c("dhsprovfrom", "dhsprovto", "riverdist")) %>%
  dplyr::rename(item1 = dhsprovfrom,
                item2 = dhsprovto)

prov.dist <- dplyr::left_join(x = prov.gcdist,
                              y = prov.roaddist, by = c("item1", "item2")) %>%
  dplyr::left_join(., y = prov.riverdist, by = c("item1", "item2")) %>%
  dplyr::rename(riverdistance = riverdist)



#..............................................................
# Combine outcome and covariate data
#..............................................................
ibDdist.prov.between$deltadistance <- purrr::pmap(ibDdist.prov.between[,c("btwn", "distcat")],
                                                  function(btwn, distcat){
  prov1 <- stringr::str_split_fixed(btwn, "/", n=2)[,1]
  prov2 <- stringr::str_split_fixed(btwn, "/", n=2)[,2]

  if(prov1 == prov2){
    deltadist <- 0
  } else {
    deltadist <- prov.dist %>%
      dplyr::filter(c(item1 == prov1 & item2 == prov2 |
                        item1 == prov2 & item2 == prov1)) %>%
      dplyr::select(distcat)
    deltadist <- unlist(unname(deltadist))
  }
  return(deltadist)
})

ibDdist.prov.between <- ibDdist.prov.between %>%
  tidyr::unnest(cols = "deltadistance")


# Scale Pred Mat and distances
scalevars <- c("prev", "precip", "temp", "elev", "crops", "netuse", "housing",
               "urban", "deltadistance")
ibDdist.prov.between <- ibDdist.prov.between %>%
  tidyr::unnest(cols = "predMat")


ibDdist.prov.between[, scalevars] <-
  apply(ibDdist.prov.between[, scalevars], 2, my.scale)


# beta for same prov or not (this is a parameter that is like a "random effect" to mop up variance for w/in prov var)
ibDdist.prov.between$provlvl <- purrr::map(ibDdist.prov.between$btwn, function(x){
  prov1 <- stringr::str_split_fixed(x, "/", n=2)[,1]
  prov2 <- stringr::str_split_fixed(x, "/", n=2)[,2]
  return(as.numeric(prov1 == prov2))
})
ibDdist.prov.between <- ibDdist.prov.between %>%
  tidyr::unnest(cols = provlvl)

#-------------------------------------------------------------------------
# CARBayes Model
#-------------------------------------------------------------------------
ibDdist.prov.between <- ibDdist.prov.between %>%
  dplyr::select(c("distcat", "F_between_mean", "provlvl", scalevars))

# group
ibDdist.prov.between <- ibDdist.prov.between %>%
  dplyr::group_by(distcat) %>%
  tidyr::nest()

#........................
# Setup here is to just
# get model form for dredges
# this model fit isn't important
#
# What is important is that we want every
# model to include the delta_distance covariate and the provlvl covariate
#.........................
options(na.action = "na.fail")
dat <- ibDdist.prov.between$data[[1]]
covar.names <- scalevars[! scalevars %in% c("deltadistance")]

mod <- lm(
  as.formula(paste("F_between_mean", "~", paste(covar.names, collapse = "+"))),
  data = dat)

# dredge
mods <- MuMIn::dredge(mod, evaluate = F,
                      m.lim = c(1,10))

# formula to string manipulation
mods <- lapply(mods, function(x){
  charform <- paste(deparse(x), collapse = "")
  ret <- str_match(charform, "F_between_mean (.*?) 1")[1,1]
  ret <- gsub("+ 1", " provlvl + deltadistance + 1", ret)
  ret <- as.formula(ret)
  return(ret)
})

# add in intercept only models
mods <- append(mods, as.formula("F_between_mean ~ provlvl + deltadistance + 1"))


#..............................................................
# Run glm model in parallel on slurm
#..............................................................
wrap_glm <- function(distcat, outcome, formula, data, burnin, n.sample){

  #..............................................................
  # fit model
  #..............................................................
  formvec <- paste(deparse(formula), collapse = "")
  betacount <- stringr::str_count(formvec, "\\+") + 1 # plus one for intercept
  prior.var.betavec <- rep(5e4, betacount) # note prior setting here

  ret <- CARBayes::S.glm(formula = formula,
                         family = "gaussian",
                         data = data,
                         burnin = burnin,
                         prior.var.beta = prior.var.betavec,
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
mod.framework.sp <- tibble(formula = mods,
                           burnin = 1e4,
                           n.sample = 1e6 + 1e4,
                           outcome = "F_between_mean"
                           )

distcat <- sort(rep(ibDdist.prov.between$distcat, length(mods)))
mod.framework.sp <- cbind.data.frame(distcat, rbind.data.frame(mod.framework.sp,
                                                                mod.framework.sp,
                                                                mod.framework.sp))

mod.framework.sp <- dplyr::left_join(mod.framework.sp, ibDdist.prov.between,
                                     by = "distcat")


# for slurm on LL
dir.create("results/carbayes_between_prov_models", recursive = T)
setwd("results/carbayes_between_prov_models/")
ntry <- 1028 # max number of nodes
sjob <- rslurm::slurm_apply(f = wrap_glm,
                            params = mod.framework.sp,
                            jobname = 'betweenMods_DICs',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = 16000,
                                                 array = sprintf("0-%d%%%d",
                                                                 ntry,
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "12:00:00"))

cat("*************************** \n Submitted Carbayes Models \n *************************** ")




################################################################################
##############             Long Final Models         ###########################
################################################################################
wrap_glm_long <- function(distcat, outcome, formula, data, burnin, n.sample){

  #..............................................................
  # fit model
  #..............................................................
  formvec <- paste(deparse(formula), collapse = "")
  betacount <- stringr::str_count(formvec, "\\+") + 1 # plus one for intercept
  prior.var.betavec <- rep(5e4, betacount) # note prior setting here

  mod <- CARBayes::S.glm(formula = formula,
                         family = "gaussian",
                         data = data,
                         burnin = burnin,
                         prior.var.beta = prior.var.betavec,
                         n.sample = n.sample)
  # note, rho by default is NULL which allows it to be estimated in the model


  #-------------------------------------------------------------------------
  # MCMC Diagnostics
  #-------------------------------------------------------------------------
  ret <- tibble::tibble(MCMC = list(mod))
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


  out <- list(
    mod = mod,
    diagnostics = ret
  )

  return(out)

}



#..............................................................
# Longer run for final models with best DIC based
# on above fits
#..............................................................
finalformulas <- list(as.formula("F_between_mean ~ provlvl + deltadistance + 1"), # gc
                      as.formula("F_between_mean ~ provlvl + deltadistance + 1"), # road
                      as.formula("F_between_mean ~ provlvl + deltadistance + 1") # river
                      )
long.mod.framework <- mod.framework.sp <- tibble(formula = finalformulas,
                                                 burnin = 1e4,
                                                 n.sample = 1e7 + 1e4,
                                                 outcome = "F_between_mean")

long.mod.framework <- cbind.data.frame(long.mod.framework, ibDdist.prov.between)

#..............................................................
# Run
#..............................................................
# we already set directory above to results/carbayes_between_prov_models
dir.create("carbayes_long_chains_BETWEEN")
mods <- purrr::pmap(long.mod.framework, wrap_glm_long)
saveRDS(mods, file = "carbayes_long_chains_BETWEEN.RDS")

