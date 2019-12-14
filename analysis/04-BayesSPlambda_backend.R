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
# Fit the lambda distributions by prov
#..............................................................

ibDdist.prov.lambda <- ibDdist.long.mtdt %>%
  dplyr::rename(riverdistance = riverdist) %>%
  dplyr::select(c("smpl1", "smpl2", "adm1name.x", "adm1name.y", "malecotf", "malecotf_gens", dplyr::ends_with("distance"))) %>%
  tidyr::gather(., key = "distcat", value = "geodist", 7:9) %>%
  dplyr::group_by(adm1name.x, distcat) %>%
  tidyr::nest()


ibDdist.prov.lambda$genbase2_fit <- purrr::map(ibDdist.prov.lambda$data, function(dat){
  dat <- dat %>%
    dplyr::filter(malecotf_gens != Inf)
  ret <- glm(malecotf_gens ~ geodist,
             data = dat,
             family = gaussian(link = "log"))
  return(ret)
})

ibDdist.prov.lambda$genbase2_slope <- purrr::map(ibDdist.prov.lambda$genbase2_fit, function(fit){
  ret <- broom::tidy(fit) %>%
    dplyr::filter(term == "geodist") %>%
    dplyr::select("estimate") %>%
    unlist(.)
  return(ret)
})

ibDdist.prov.lambda <- ibDdist.prov.lambda %>%
  tidyr::unnest(cols = "genbase2_slope")


#..............................................................
# Get within IBD for the Prov
#..............................................................
ibDdist.prov.lambda$withinprovIBD <- purrr::map(ibDdist.prov.lambda$data, function(x){
  ret <- x %>%
    dplyr::filter(geodist == 0) %>%
    dplyr::filter(malecotf_gens != Inf) %>%
    dplyr::summarise(
      meangens = mean(malecotf_gens)
    )
  return(as.numeric(ret))
})

ibDdist.prov.lambda <- ibDdist.prov.lambda %>%
  tidyr::unnest(cols = c("withinprovIBD"))


#..............................................................
# Get mean IBD for the Prov
#..............................................................
ibDdist.prov.lambda$meanprovIBD <- purrr::map(ibDdist.prov.lambda$data, function(x){
  ret <- x %>%
    dplyr::filter(malecotf_gens != Inf) %>%
    dplyr::summarise(
      meangens = mean(malecotf_gens)
    )
  return(as.numeric(ret))
})

ibDdist.prov.lambda <- ibDdist.prov.lambda %>%
  tidyr::unnest(cols = c("meanprovIBD"))


#..............................................................
# Import Covariates for Province
#..............................................................
provCovar <- readRDS("data/derived_data/covar_rasterstack_provlocations_scaled.RDS")

#..............................................................
# Combine outcome and covariate data
#..............................................................
# housekeeping
riskvar <- colnames(provCovar)[!colnames(provCovar) %in% "adm1name"]
ibDdist.prov.lambda <- ibDdist.prov.lambda %>%
  dplyr::rename(adm1name = adm1name.x) %>%
  dplyr::select(c("adm1name", "distcat", "genbase2_slope", "withinprovIBD", "meanprovIBD"))


# model framework
mod.IBD.provCovar <- dplyr::left_join(ibDdist.prov.lambda, provCovar,
                                       by = "adm1name")


mod.IBD.provCovar <- mod.IBD.provCovar %>%
  dplyr::select("adm1name", "distcat", riskvar, "genbase2_slope", "withinprovIBD", "meanprovIBD")

mod.IBD.provCovar.nest <- mod.IBD.provCovar %>%
  tidyr::gather(., key = "outcome", value = "ibd", 12:14) %>%
  dplyr::group_by(distcat, outcome) %>%
  tidyr::nest()


#-------------------------------------------------------------------------
# Conditional Autoregressive Spatial Model
#-------------------------------------------------------------------------
#......................
# Make Adjacency Matrix for Pv
#......................
W.nb <- spdep::poly2nb(sf::as_Spatial(DRCprov), row.names = DRCprov$adm1name)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)

#..............................................................
# Make Model Framework
#..............................................................
prov.covar.names <- names(provCovar)[names(provCovar) != "adm1name"]

#........................
# Setup here is to just
# get model form for dredges
# this model fit isn't important
#.........................
options(na.action = "na.fail")
dat <- mod.IBD.provCovar.nest$data[[1]]
mod <- lm(
  as.formula(paste("ibd", "~", paste(prov.covar.names, collapse = "+"))),
  data = dat)

# dredge
mods <- MuMIn::dredge(mod, evaluate = F,
                      m.lim = c(1,10))

# formula to string manipulation
mods <- lapply(mods, function(x){
  charform <- paste(deparse(x), collapse = "")
  ret <- str_match(charform, "ibd (.*?) 1")[1,1]
  ret <- as.formula(ret)
  return(ret)
})

# add in intercept only models
mods <- append(mods, as.formula("ibd ~ 1"))


#..............................................................
# Run Leroux model in parallel on slurm
#..............................................................
wrap_S.CARleroux <- function(distcat, formula, W, data, burnin, n.sample){

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
                           W = list(W))

# rep this out three times for levels of spatial data
mod.framework.sp <- lapply(1:3, function(x) return(mod.framework.sp)) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(distcat = c( rep("gcdistance", times = nrow(mod.framework.sp)),
                             rep("roaddistance", times = nrow(mod.framework.sp)),
                             rep("riverdistance", times = nrow(mod.framework.sp))
                          )
  )

mod.IBD.provCovar.nest.framework.sp <- dplyr::left_join(mod.IBD.provCovar.nest,
                                                        mod.framework.sp,
                                                        by = "distcat")


# for slurm on LL
dir.create("results/carbayes_sp_dics", recursive = T)
setwd("results/carbayes_sp_dics/")
ntry <- 1028 # max number of nodes
sjob <- rslurm::slurm_apply(f = wrap_S.CARleroux,
                            params = mod.IBD.provCovar.nest,
                            jobname = 'CARleroux_DICs',
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

