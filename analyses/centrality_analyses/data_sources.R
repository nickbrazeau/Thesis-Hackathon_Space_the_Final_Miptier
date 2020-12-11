####################################################################################
## Purpose: Pull together data sources for python import
##
##
####################################################################################

library(tidyverse)
dir.create("data/python_import")

#............................................................
# geodistances
#...........................................................
# symmetric matrix 492 x 492 for the 492 "villages" (clusters) --> in long format
# Not all clusters have malaria genetic samples
geodists <- readRDS("data/distance_data/distancematrix_bycluster.rds")
readr::write_csv(geodists, "data/python_import/distancematrix_bycluster.csv")


#............................................................
# gradient descent results
#...........................................................
finbd <- readRDS("results/clust_inbd_results/all_smpls/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS")
finbd_clust <- finbd %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = "inbreed_ests") %>%
  dplyr::filter(param != "m") %>%
  dplyr::rename(hv001 = param,
                Finbd = est)

readr::write_csv(finbd_clust, "data/python_import/cluster_inbreedinf_coefficient_results_351x.csv")


#............................................................
# cluster covars
# ubranicity and pf-incidence
#...........................................................
clust.covars <- readRDS("data/derived_data/covar_rasterstack_samplinglocations_raw.RDS")
clust.covars <- clust.covars %>% dplyr::select(c("hv001", "incidence", "urban"))
# long-lat for clusters
ge <- readRDS("data/derived_data/spacemips_GE.rds") %>%
  dplyr::select(-c("geometry"))
clust.covars <- dplyr::inner_join(ge, clust.covars, by = "hv001")
readr::write_csv(clust.covars, "data/python_import/clust_covariates_351x.csv")

#......................
# pop density
#......................
popdensdf <- readRDS("data/derived_data/popdens_samplinglocations_raw.RDS") %>%
  dplyr::select(-c("worldpop_pop_dens", "geometry"))
sf::st_geometry(popdensdf) <- NULL
readr::write_csv(popdensdf, "data/python_import/pop_density.csv")


#............................................................
# highly related pairs
#...........................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
ibD.meiotics <- ibD %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf")) %>%
  dplyr::filter(malecotf >= 0.5)
# output csv files for importing into python
readr::write_csv(ibD.meiotics, "data/python_import/meioticse_ibD_long_mtdt.csv")

#......................
# cluster counts for hightly related pairs
#......................
ge <- ge %>%
  dplyr::filter(hv001 %in% clust.covars$hv001)

withins <- ibD.meiotics %>%
  dplyr::filter(hv001.x == hv001.y)
within_clsts <- table(withins$hv001.x)
within_clsts <- tibble::tibble(
  hv001 = names(within_clsts),
  counts = unname(unlist(within_clsts)))


# betweens
betweens <- ibD.meiotics %>%
  dplyr::filter(hv001.x != hv001.y)

betweens_clsts <- c(betweens$hv001.x, betweens$hv001.y)
betweens_clsts <- table(betweens_clsts)
betweens_clsts <- tibble::tibble(
  hv001 = names(betweens_clsts),
  counts = unname(unlist(betweens_clsts)))


# output csv file fo rbetween for importing into python
readr::write_csv(betweens_clsts, "data/python_import/between_highly_relat_cluster_counts.csv")
readr::write_csv(within_clsts, "data/python_import/within_highly_relat_cluster_counts.csv")
