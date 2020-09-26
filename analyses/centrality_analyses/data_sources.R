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

# asymmetric matrix 38 x 38 for migration rate flows between
# voroni tesselated "territories" or provinces --> long format

# migrateflows <- readRDS("data/distance_data/voroni_migration_flows_fromworldpop.RDS") %>%
migrateflows <- readRDS("data/distance_data/vr_nodepairs_migrate_disance.rds") %>%
  dplyr::select(c("NODEI", "NODEJ", "PrdMIG", "PrdMIG_scaled")) %>%
  dplyr::filter(!duplicated(.))
readr::write_csv(migrateflows, "data/python_import/voroni_migration_flows_fromworldpop.csv")


#............................................................
# gradient descent results
#...........................................................
finbd <- readRDS("results/clust_inbd_results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS")
finbd_clust <- finbd %>%
  dplyr::filter(spacetype != "migrate") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = "inbreed_ests") %>%
  dplyr::filter(param != "m") %>%
  dplyr::rename(hv001 = param,
                Finbd = est)

readr::write_csv(finbd_clust, "data/python_import/cluster_inbreedinf_coefficient_results_351x.csv")


# provinces
finbd_prov <- finbd %>%
  dplyr::filter(spacetype == "migrate") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = "inbreed_ests") %>%
  dplyr::filter(param != "m") %>%
  dplyr::rename(hv001 = param,
                Finbd = est)
readr::write_csv(finbd_prov, "data/python_import/voroniterritories_inbreedinf_coefficient_results_38x.csv")



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

#...........................................................
# voroni migration territories covars ("province")
# ubranicity and pf-incidence
#...........................................................
prov.covars <- readRDS("data/derived_data/covar_rasterstack_provterritories_raw.RDS")
prov.covars <- prov.covars %>% dplyr::select(c("IPUMSID", "incidence", "urban"))
vrdf <- readRDS("data/distance_data/voroni_base.RDS")
prov.covars <- dplyr::left_join(vrdf, prov.covars, by = "IPUMSID") %>%
  dplyr::rename(geometry = x) %>%
  dplyr::mutate(centroid = sf::st_centroid(geometry),
                centroidlong = sf::st_coordinates(centroid)[,1],
                centroidlat = sf::st_coordinates(centroid)[,2]) %>%
  dplyr::select(-c("centroid"))
readr::write_csv(prov.covars, "data/python_import/voroniterritories_covariates_38x.csv")



#............................................................
# highly related pairs
#...........................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
ibD.meiotics <- ibD %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf")) %>%
  dplyr::filter(malecotf >= 0.5)

# output csv files for importing into python 


readr::write_csv(ibD, "data/python_import/raw_pairwise_ibD_long_mtdt.csv")

