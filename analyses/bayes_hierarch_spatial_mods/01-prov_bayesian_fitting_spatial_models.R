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
####  Areal-Based Modeling ####
#......................................................................................................
library(CARBayes)
#......................
# bring in data
#......................
# covariate data
prov_covars <- readRDS("data/derived_data/covar_rasterstack_provterritories_scaled.RDS") %>%
  dplyr::mutate(param = as.character(IPUMSID))

# load pretty map aesthetics
load("data/map_bases/space_mips_maps_bases.rda")

# spatial data
vrdf <- readRDS("data/distance_data/voroni_base.RDS") %>%
  dplyr::rename(param = IPUMSID) %>%
  dplyr::mutate(param = as.character(param)) %>%
  rename(geometry = x)

# inbreeding data
prov_inbd <- readRDS("results/clust_inbd_results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS") %>%
  dplyr::filter(spacetype == "migrate") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = inbreed_ests) %>%
  dplyr::mutate(logitest = logit(est)) %>%
  dplyr::ungroup(.)

# bring together
prov_spatdat <- dplyr::left_join(prov_inbd, prov_covars, by = "param") %>%
  dplyr::left_join(., vrdf, by = "param") %>%
  dplyr::filter(param != "m")
prov_spatdat <- sf::st_as_sf(prov_spatdat)
sf::st_crs(prov_spatdat) <-  "+proj=longlat +datum=WGS84 +no_defs"


#......................
# Make Adjacency Matrix for Pv
#......................
W.nb <- spdep::poly2nb(prov_spatdat, row.names = prov_spatdat$param)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)

#......................
# Make Model Framework
#......................
prov_mod_framework <- tibble(name = c("CAR_intercept",
                                      "interaction"
),
formula = c("logitest ~ 1",
            "logitest ~ incidence + urban + incidence*urban"
),
burnin = 1e4,
n.sample = 1e4 + 1e4)

# add needed bits to model
prov_mod_framework$data <- lapply(1:nrow(prov_mod_framework), function(x) return(prov_spatdat))
prov_mod_framework$W <- lapply(1:nrow(prov_mod_framework), function(x) return(W))

# Make a wrapper for leroux way
wrap_S.CARleroux <- function(name, formula, W, data, burnin, n.sample){

  formvec <- paste(formula, collapse = "")
  betacount <- stringr::str_count(formvec, "\\+") + 2 # need intercept and betas
  betacount <- ifelse(grepl("1", formvec), 1, betacount) # corner case of just intercept
  prior.var.betavec <- rep(5e4, betacount) # note prior setting here

  # rho here is NULL by default and esimtated in the model
  ret <- CARBayes::S.CARleroux(formula = as.formula(formula),
                               family = "gaussian",
                               thin = 100,
                               W = W,
                               data = data,
                               burnin = burnin,
                               prior.var.beta = prior.var.betavec,
                               n.sample = n.sample)


  return(ret)
}


prov_mod_framework$MCMC <- purrr::pmap(prov_mod_framework, wrap_S.CARleroux)

#......................
# get posteriors
#......................
prov_mod_framework$mcmc.modsum <- purrr::map(prov_mod_framework$MCMC, print) # note, print is overloaded here
prov_mod_framework$summresults <- purrr::map(prov_mod_framework$mcmc.modsum, "summary.results")
# get samples for HPD
prov_mod_framework <- prov_mod_framework %>%
  dplyr::mutate(samples = purrr::map(MCMC, "samples"))
# get HPD
prov_mod_framework$HPD <- purrr::map(prov_mod_framework$samples, HDInterval::hdi, credMass = 0.95)
# subset HPD to relevant params
prov_mod_framework$HPD.summ <- purrr::map(prov_mod_framework$HPD,
                                          function(x){
                                            ret <- x[ names(x) %in% c("beta", "rho", "nu2", "tau2") ]
                                            ret <- lapply(ret, function(y){
                                              ret <- t(y)
                                              colnames(ret) <- c("HPD 0.025", "HPD 0.975")
                                              return(ret) }) %>%
                                              do.call("rbind.data.frame", .)
                                          })

prov_mod_framework$summresultsdf <- purrr::pmap(prov_mod_framework[,c("summresults", "HPD.summ")],
                                                function(summresults, HPD.summ){
                                                  pars <- rownames(summresults)
                                                  ret <- cbind.data.frame(pars = pars,
                                                                          cbind( summresults, HPD.summ))
                                                  return(ret)
                                                })

# Means for comparison to median
prov_mod_framework$posteriormeans <- purrr::pmap(prov_mod_framework[,c("samples", "summresults")],
                                                 function(samples, summresults){
                                                   betas <- apply(samples$beta, 2, function(x){mean(expit(x))})
                                                   names(betas) <- row.names(summresults)[!row.names(summresults) %in% c("tau2", "nu2", "rho")]
                                                   tau2 <-  apply(samples$tau2, 2, function(x){mean(expit(x))})
                                                   names(tau2) <- "tau2"
                                                   rho <- apply(samples$rho, 2, function(x){mean(expit(x))})
                                                   names(rho) <- "rho"

                                                   params <- matrix(NA, nrow = 1, ncol = length(c(betas, tau2, rho)))
                                                   params <- c(betas, tau2, rho)
                                                   params <- data.frame(covar = names(params), postmean = params)
                                                   return(params)
                                                 })

# also want to get phis here for each province
prov_mod_framework$posteriorphis <- purrr::map(prov_mod_framework$samples,
                                               function(smpl){

                                                 phimean <- apply(smpl$phi, 2, function(x){mean(expit(x))})
                                                 phimedian <- apply(smpl$phi, 2, function(x){median(expit(x))})
                                                 phiperc025 <- apply(smpl$phi, 2, function(x){quantile(expit(x), 0.025)})
                                                 phiperc975 <- apply(smpl$phi, 2, function(x){quantile(expit(x), 0.975)})

                                                 phihpd <- HDInterval::hdi(expit(smpl$phi), credMass = 0.95)

                                                 phineff <- apply(smpl$phi, 2, coda::effectiveSize)
                                                 phis <- list(phimean, phimedian, phiperc025, phiperc975, phihpd, phineff)
                                                 phis <- do.call("rbind.data.frame", phis)
                                                 # bind together
                                                 phis <- cbind.data.frame(statistic = c("Mean", "Median", "Perc0.025", "Perc97.5",
                                                                                        "HPD 0.025", "HPD 0.975", "Neff"),
                                                                          phis)
                                                 # put names in
                                                 colnames(phis)[2:ncol(phis)] <- prov_spatdat$param
                                                 return(phis)
                                               })


prov_mod_framework_res <- prov_mod_framework %>%
  dplyr::select(c("name", "formula", "summresultsdf")) %>%
  tidyr::unnest(cols = "summresultsdf")


#......................
# DICs
#......................
prov_mod_framework$modfit <- purrr::map(prov_mod_framework$mcmc.modsum, "modelfit")

# note in order to be consistent, with the PrevMap approach,
# I want Gelman's DIC: DICg = mu + sigma
# can get this from the log likelihood provided by CarBayes.

# CarBayes nicely outputs the Fitted Values, from which I can calculate a log-likelihood
dimnames(prov_mod_framework$W[[1]])[[1]]
# prov_mod_framework$LL <- purrr::pmap(prov_mod_framework[, c("data", "trials", "MCMC")],
#                                      function(data, trials, MCMC){
#                                        Yi <- data$plsmdn
#                                        ni <- trials
#                                        Ti <- MCMC$samples$fitted
#
#                                        LL <- c()
#                                        for(s in 1:nrow(Ti)){ # for every sim
#                                          LL.iter <- 0
#                                          for(i in 1:ncol(Ti)){ # for every spatial unit
#                                            LL.iter <- dbinom(x = Yi[i], size = ni[i], prob = Ti[s, i]/ni[i], log = T) + LL.iter
#                                            LL.iter
#                                          }
#                                          LL <- append(LL.iter, LL)
#
#                                        }
#
#                                        LL <- return(LL)
#                                      })


prov_mod_framework$DICg <- purrr::map(prov_mod_framework$LL, function(x){
  mu = mean(x)
  sigma = var(x)/2
  DICg = mu + sigma
  return(DICg)
})


prov_mod_framework$DICs <- purrr::map(prov_mod_framework$modfit, "DIC")

#............................................................
# predictions
#...........................................................
provpred <- prov_mod_framework
provpred <- provpred %>%
  dplyr::mutate(fittedvalues = purrr::map(MCMC, "samples"),
                fittedvalues = purrr::map(fittedvalues, "fitted"),
                provmodel.posteriormeans = purrr::map(fittedvalues,
                                                    function(x){apply(x, 2, function(y){mean(expit(y))})}),
                provmodel.posteriormeans = purrr::map(provmodel.posteriormeans, function(x){cbind.data.frame(param = as.character(provpred$data[[1]]$IPUMSID),
                                                                                                             posteriormeans = x)}))
provpred <- provpred %>%
  dplyr::select(c("name", "provmodel.posteriormeans"))


#......................
# plot out
#......................
# intercept model
intercept_provmodelplotObjdf <- dplyr::left_join(vrdf,
                                                 provpred$provmodel.posteriormeans[provpred$name == "CAR_intercept"][[1]])
intercept_provmodelplotObj <- ggplot() +
  geom_sf(data = intercept_provmodelplotObjdf, aes(fill = posteriormeans)) +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1) +
  prettybasemap_nodrc_nonorth_dark +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0),
    legend.text = element_text(face = "bold", size = 11),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))



# interaction model
interaction_provmodelplotObjdf <- dplyr::left_join(vrdf,
                                                 provpred$provmodel.posteriormeans[provpred$name == "interaction"][[1]])
interaction_provmodelplotObj <- ggplot() +
  geom_sf(data = interaction_provmodelplotObjdf, aes(fill = posteriormeans)) +
  scale_fill_viridis_c("Inbreeding", option="plasma", direction = 1) +
  prettybasemap_nodrc_nonorth_dark +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12, vjust = 0.5, hjust = 0),
    legend.text = element_text(face = "bold", size = 11),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05),"cm"))


# save
saveRDS(prov_mod_framework, file = "results/clust_inbd_results/final_prov_maps/bayesian_posterior_prov_inbd_vrdf_terrmap_fits.RDS")
saveRDS(provpred, file = "results/clust_inbd_results/final_prov_maps/bayesian_posterior_prov_inbd_vrdf_terrmap_predictions.RDS")

saveRDS(intercept_provmodelplotObj, file = "results/clust_inbd_results/final_prov_maps/bayesian_posterior_prov_inbd_intercept_plot.RDS")
saveRDS(interaction_provmodelplotObj, file = "results/clust_inbd_results/final_prov_maps/bayesian_posterior_prov_inbd_interaction_plot.RDS")




