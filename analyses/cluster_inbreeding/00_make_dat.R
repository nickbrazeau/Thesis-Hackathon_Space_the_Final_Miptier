#########################################################################
# Purpose: Wrangle data for cluster inbreeding grad descent function
#
# Date: April 02 2020
#########################################################################
source("R/pairwise_helpers.R")
library(tidyverse)

#......................
# wrangle and take to long, diagonals to 0 for distances
#......................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster) %>%
  dplyr::mutate(hv001.x = as.numeric(as.character(hv001.x)),
                hv001.y = as.numeric(as.character(hv001.y))) # factor liftover

gcdist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "gcdistance")) %>%
  dplyr::mutate(gcdistance = ifelse(hv001.x == hv001.y, 0, gcdistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


roaddist_gens <- ibD %>%
  dplyr::left_join(., distancematrix.cluster, by = c("hv001.x", "hv001.y")) %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf", "roaddistance")) %>%
  dplyr::mutate(roaddistance = ifelse(hv001.x == hv001.y, 0, roaddistance)) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))

#......................
# wrangle migration provinces demes
#......................
# note can't use expand_distance_matrix here bc assumes symmetrical distance matrix -- had to deal w/ upstream
node_pairs <- readRDS(file = "data/distance_data/vr_nodepairs_migrate_disance.rds")
# note distancematrix expansion doesn't include clut = clust, so just doing the same self is 0 distance as above, so need to have sames
extra_nodes <- readRDS("data/distance_data/vr_memberships.rds") %>%
  dplyr::select(-c("x")) %>%
  dplyr::rename(hv001.x = hv001,
                NODEI = IPUMSID) %>%
  dplyr::mutate(hv001.y = hv001.x,
                NODEJ = NODEI,
                ISO = NA, # extra bits for binding above
                LONFR = NA,
                LATFR = NA,
                LONTO = NA,
                LATTO = NA,
                PrdMIG = 0,
                PrdMIG_scaled = 0)
# bring together
node_pairs <- dplyr::bind_rows(node_pairs, extra_nodes)


# need to expand IBD here b/c of asymetry again
ibD <- ibD %>%
  dplyr::select(c("smpl1", "smpl2", "hv001.x", "hv001.y", "malecotf"))
ibD.copy <- ibD
colnames(ibD.copy) <- c("smpl2", "smpl1", "hv001.y", "hv001.x", "malecotf")

# NB, duplication filter in the Cpp inbreeding coeff calculator will
# first allow this to be replicated (for usual symetric distances)
# and then filter the duplicates
ibD.long <- dplyr::bind_rows(ibD, ibD.copy)
migrate_gens <- ibD.long %>%
  dplyr::left_join(., node_pairs, by = c("hv001.x", "hv001.y")) %>%
  dplyr::mutate(PrdMIG = ifelse(NODEI == NODEJ, 0, PrdMIG)) %>%
  dplyr::select(c("smpl1", "smpl2", "NODEI", "NODEJ", "malecotf", "PrdMIG")) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))


#......................
# save out
#......................
dir.create("data/derived_data/clst_inbreeding_dat")
saveRDS(as.data.frame(gcdist_gens),
        file = "data/derived_data/clst_inbreeding_dat/gcdist_gens.RDS")
saveRDS(as.data.frame(roaddist_gens),
        file = "data/derived_data/clst_inbreeding_dat/roaddist_gens.RDS")
saveRDS(as.data.frame(migrate_gens),
        file = "data/derived_data/clst_inbreeding_dat/migrate_gens.RDS")


