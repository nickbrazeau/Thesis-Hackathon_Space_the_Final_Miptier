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
finbd <- readRDS("results/clust_inbd_results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS")
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
ge <- readRDS("data/derived_data/spacemips_GE.rds")
clust.covars <- dplyr::left_join(ge, clust.covars, by = "hv001")
readr::write_csv(clust.covars, "data/python_import/clust_covariates_351x.csv")

#............................................................
# highly related pairs
#...........................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
ibD.meiotics <- ibD %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf")) %>%
  dplyr::filter(malecotf >= 0.5)

# output csv files for importing into python


readr::write_csv(ibD, "data/python_import/raw_pairwise_ibD_long_mtdt.csv")

