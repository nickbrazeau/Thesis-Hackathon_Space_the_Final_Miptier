#..............................................................
# Purpose of this script is to fit the various models
# from CarBayes for our IBD Spatial Covariates.
# We use this "backend" script to interact with a server
#
# Y_k ~ B^T_X_{pred} + O_k + \phi_k
#..............................................................
library(tidyverse)
library(CARBayes)
library(raster)
source("R/basics.R")
source("R/pairwise_helpers.R")
set.seed(48)

# drc prov for plot
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))


#....................................................................................
# Import Genetic Data
#....................................................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
# log2 IBD measure
ibD <- ibD %>%
  dplyr::mutate(malecotf_gens = -log2(malecotf),
                malecotf_gens_inv = 1/malecotf_gens)

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


################################################################################
##############             Diagnostic Models         ###########################
################################################################################
#..............................................................
# Get within IBD for the Prov
#..............................................................
# here we want the long format since we want to be able to group_by
# adm1name.x so need every pairwise to be there
ibDdist.prov.within <- ibDdist.long %>%
  dplyr::rename(riverdistance = riverdist) %>%
  dplyr::select(c("smpl1", "smpl2", "adm1name.x", "adm1name.y", "malecotf", "malecotf_gens", dplyr::ends_with("distance"))) %>%
  tidyr::gather(., key = "distcat", value = "geodist", 7:9) %>%
  dplyr::group_by(adm1name.x, distcat) %>%
  tidyr::nest()

ibDdist.prov.within$withinprovIBD <- purrr::map(ibDdist.prov.within$data, function(x){
  ret <- x %>%
    dplyr::filter(geodist == 0) %>%
    dplyr::summarise(
      meangens = mean(malecotf)
    )
  return(as.numeric(ret))
})

ibDdist.prov.within <- ibDdist.prov.within %>%
  tidyr::unnest(cols = c("withinprovIBD"))



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
    actuse = my.scale(logit(actuse, tol = 0.1)),
    netuse = my.scale(logit(netuse, tol = 0.1)),
    housing = my.scale(logit(housing, tol = 0.1))
  )


#..............................................................
# Combine outcome and covariate data
#..............................................................
# housekeeping
riskvar <- colnames(provCovar)[!colnames(provCovar) %in% "adm1name"]
ibDdist.prov.within <- ibDdist.prov.within %>%
  dplyr::rename(adm1name = adm1name.x) %>%
  dplyr::select(c("adm1name", "distcat", "withinprovIBD"))


# model framework
mod.IBD.provCovar <- dplyr::left_join(ibDdist.prov.within, provCovar,
                                       by = "adm1name") %>%
  dplyr::select("adm1name", "distcat", riskvar, "withinprovIBD")

mod.IBD.provCovar.nest <- mod.IBD.provCovar %>%
  dplyr::group_by(distcat) %>%
  tidyr::nest()


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
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
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
W.gcdist <- make_symm_mat(W.gcdist)
scale <- mean(W.gcdist[lower.tri(W.gcdist, diag = T)])
W.gcdist <- dexp(W.gcdist/scale) # scale to help with variance and change rate


# road
prov.roaddist <- readRDS("data/distance_data/prov_road_distmeters_long.rds")
W.roaddist <- prov.roaddist %>%
  tidyr::spread(., key = "item2", value = "roaddistance") %>%
  dplyr::select(-c("item1"))
# now make it a symm matrix
W.roaddist <- make_symm_mat(W.roaddist)
scale <- mean(W.roaddist[lower.tri(W.roaddist, diag = T)])
W.roaddist <- dexp(W.roaddist/scale) # scale to help with variance and change rate

# river
prov.riverdist <- readRDS("data/distance_data/river_distance_forprovinces.rds")
W.riverdist <- prov.riverdist %>%
  dplyr::select(c("dhsprovfrom", "dhsprovto", "riverdist")) %>%
  dplyr::rename(item1 = dhsprovfrom,
                item2 = dhsprovto) %>%
  tidyr::spread(., key = "item2", value = "riverdist") %>%
  dplyr::select(-c("item1"))
# now make it a symm matrix
W.riverdist <- make_symm_mat(W.riverdist)
scale <- mean(W.riverdist[lower.tri(W.riverdist, diag = T)])
W.riverdist <- dexp(W.riverdist/scale) # scale to help with variance and change rate

#..............................................................
# Make Model Framework
#..............................................................
prov.covar.names <- names(provCovar)[names(provCovar) != "adm1name"]

#........................
# Setup here is to just
# get model form for dredges
# this model fit isn't important
#.........................
covars.eval <- c("prev", "precip", "temp", "elev", "crops", "netuse", "housing", "urban")
options(na.action = "na.fail")
dat <- mod.IBD.provCovar.nest$data[[1]]
mod <- lm(
  as.formula(paste("withinprovIBD", "~", paste(covars.eval, collapse = "+"))),
  data = dat)

# dredge
mods <- MuMIn::dredge(mod, evaluate = F,
                      m.lim = c(1,10))

# formula to string manipulation
mods <- lapply(mods, function(x){
  charform <- paste(deparse(x), collapse = "")
  ret <- str_match(charform, "withinprovIBD (.*?) 1")[1,1]
  ret <- as.formula(ret)
  return(ret)
})

# add in intercept only models
mods <- append(mods, as.formula("withinprovIBD ~ 1"))


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
mod.framework.sp <- tibble(formula = mods,
                           burnin = 1e4,
                           n.sample = 1e6 + 1e4,
                           outcome = "withinprovIBD")

W.matrices <- tibble::tibble(
  distcat = c("gcdistance", "roaddistance", "riverdistance"),
  W = list(W.gcdist, W.roaddist, W.riverdist))

W.matrices <- lapply(1:nrow(mod.framework.sp), function(x) return(W.matrices)) %>%
  dplyr::bind_rows(.) %>%
  dplyr::arrange(., distcat)

mod.framework.sp <- cbind.data.frame(
  rbind.data.frame(mod.framework.sp, mod.framework.sp, mod.framework.sp),
  W.matrices
)

mod.framework.sp <- dplyr::left_join(mod.framework.sp, mod.IBD.provCovar.nest,
                                     by = "distcat")


# for slurm on LL
scrdir <- "/pine/scr/n/f/nfb/Projects/Space_the_Final_Miptier/results/carbayes_within_prov_models"
dir.create(scrdir, recursive = T)
setwd(scrdir)
ntry <- nrow(mod.framework.sp) # max number of nodes
sjob <- rslurm::slurm_apply(f = wrap_S.CARleroux,
                            params = mod.framework.sp,
                            jobname = 'CARleroux_DICs_within',
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
wrap_S.CARleroux_long <- function(distcat, outcome, formula, W, data, burnin, n.sample){

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
# Longer run for final models with best DIC based
# on above fits
#..............................................................
finalformulas <- list(as.formula("withinprovIBD ~ 1"), # gc
                      as.formula("withinprovIBD ~ 1"), # road
                      as.formula("withinprovIBD ~ 1") # river
)
long.mod.framework <- mod.framework.sp <- tibble(formula = finalformulas,
                                                 burnin = 1e4,
                                                 n.sample = 1e7 + 1e4,
                                                 outcome = "withinprovIBD",
                                                 W = list(W.gcdist, W.roaddist, W.riverdist))

long.mod.framework <- cbind.data.frame(long.mod.framework, mod.IBD.provCovar.nest)

#..............................................................
# Run
#..............................................................
# we already set directory above to results/carbayes_within_prov_models
mods <- purrr::pmap(long.mod.framework, wrap_S.CARleroux_long)
saveRDS(mods, file = "carbayes_long_chains_WITHIN.RDS")





