source("R/pairwise_helpers.R")
#..............................................................
# meta
#..............................................................
drcsmpls <- readRDS("data/distance_data/drcsmpls_foruse.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("id", "hv001")) %>%
  dplyr::rename(name = id)

mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(name = id) %>%
  dplyr::select(c("name", "country", "hv001", "adm1name", "longnum", "latnum")) %>%
  dplyr::filter(name %in% drcsmpls$name)

#...................
# Genetic Data
#...................
# IBD
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD_polarized_biallelic_processed.long.rds")
# merge in cluster informations
colnames(drcsmpls)[1] <- "smpl1"
ibD <- dplyr::left_join(ibD, drcsmpls, by = "smpl1")

colnames(drcsmpls)[1] <- "smpl2"
ibD <- dplyr::left_join(ibD, drcsmpls, by = "smpl2")
# rename for recursive function
ibD <- ibD %>%
  dplyr::rename(relatedness = malecotf)


#...................
# Distance Data
#...................
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- distancematrix.cluster %>%
  dplyr::filter(hv001.x %in% unique(drcsmpls$hv001)) %>%
  dplyr::filter(hv001.y %in% unique(drcsmpls$hv001))
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster)

gcdistmat <- distancematrix.cluster %>%
  dplyr::select(c("hv001.x", "hv001.y", "gcdistance")) %>%
  dplyr::rename(distance = gcdistance) %>%
  dplyr::mutate(distance = distance/1e3)

#..............................................................
# get data
#..............................................................
distdata <- dplyr::left_join(ibD, gcdistmat, by = c("hv001.x", "hv001.y"))



#..............................................................
# downsize setup
#..............................................................
keep <- c(1:5)
smpls.keep <- drcsmpls %>%
  dplyr::filter(hv001 %in% keep) %>%
  dplyr::select(c("smpl2")) %>%
  unlist(.)


#..............................................................
# downsize
#..............................................................
distdata <- distdata %>%
  dplyr::filter(c(smpl1 %in% smpls.keep & smpl2 %in% smpls.keep ))

distdata <- distdata %>%
  dplyr::mutate(distance = ifelse(hv001.x == hv001.y, 0, distance))
colnames(distdata) <- c("smpl1", "smpl2", "relatedness", "K1", "K2", "distance")

