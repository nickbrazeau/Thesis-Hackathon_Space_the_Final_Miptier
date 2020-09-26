####################################################################################
## Purpose: Pull together data sources
##
##
####################################################################################
library(dplyr) 

#............................................................
# geodistances
#...........................................................
# symmetric matrix 492 x 492 for the 492 "villages" (clusters) --> in long format
# Not all clusters have malaria genetic samples
geodists <- readRDS("data/distance_data/distancematrix_bycluster.rds")

# asymmetric matrix 38 x 38 for migration rate flows between
# voroni tesselated "territories" or provinces --> long format

# migrateflows <- readRDS("data/distance_data/voroni_migration_flows_fromworldpop.RDS") %>%
migrateflows <- readRDS("data/distance_data/vr_nodepairs_migrate_disance.rds") %>%
  dplyr::select(c("NODEI", "NODEJ", "PrdMIG", "PrdMIG_scaled")) %>%
  dplyr::filter(!duplicated(.))

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

finbd_prov <- finbd %>%
  dplyr::filter(spacetype == "migrate") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = "inbreed_ests") %>%
  dplyr::filter(param != "m") %>%
  dplyr::rename(hv001 = param,
                Finbd = est)



#............................................................
# cluster covars
# ubranicity and pf-incidence
#...........................................................
clust.covars <- readRDS("data/derived_data/covar_rasterstack_samplinglocations_raw.RDS")
clust.covars <- clust.covars %>% dplyr::select(c("hv001", "incidence", "urban")) %>%
  dplyr::mutate(hv001 = as.character(hv001)) # num to char -- ok
# long-lat for clusters
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))
# DRC prov base
ge <- readRDS("data/derived_data/spacemips_GE.rds")

#...........................................................
# voroni migration territories covars ("province")
# ubranicity and pf-incidence
#...........................................................
prov.covars <- readRDS("data/derived_data/covar_rasterstack_provterritories_raw.RDS")
prov.covars <- prov.covars %>% dplyr::select(c("IPUMSID", "incidence", "urban"))


#............................................................
# highly related pairs
#...........................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
ibD.meiotics <- ibD %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf")) %>%
  dplyr::filter(malecotf >= 0.5)

# output csv files for importing into python 


# sanity
