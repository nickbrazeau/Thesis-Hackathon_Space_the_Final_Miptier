#..............................................................
# Purpose of this script is to calculate distance
# likelihood
#..............................................................
library(tidyverse)
source("R/pairwise_helpers.R")
source("R/recursive_likelihood.R")

#..............................................................
# Model Setup
#..............................................................
mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("name", "country", "hv001", "adm1name", "longnum", "latnum"))

#...................
# Genetic Data
#...................
# IBD
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
# rename for recursive function
ibD <- ibD %>%
  dplyr::rename(relatedness = malecotf) %>%
  dplyr::select(c("smpl1", "smpl2", "relatedness", "hv001.x", "hv001.y"))

# ibs
ibS <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS.long.mtdt.rds")
# rename for recursive function
ibS <- ibS %>%
  dplyr::rename(relatedness = hammings) %>%
  dplyr::select(c("smpl1", "smpl2", "relatedness", "hv001.x", "hv001.y"))

#...................
# Distance Data
#...................
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.cluster <- distancematrix.cluster %>%
  dplyr::filter(hv001.x %in% unique(mtdt$hv001)) %>%
  dplyr::filter(hv001.y %in% unique(mtdt$hv001))
distancematrix.cluster <- expand_distance_matrix(distancematrix.cluster)

gcdistmat <- distancematrix.cluster %>%
  dplyr::select(c("hv001.x", "hv001.y", "gcdistance")) %>%
  dplyr::rename(distance = gcdistance) %>%
  dplyr::mutate(distance = distance/1e3)

roaddistmat <- distancematrix.cluster %>%
  dplyr::select(c("hv001.x", "hv001.y", "roaddistance")) %>%
  dplyr::rename(distance = roaddistance) %>%
  dplyr::mutate(distance = distance/1e3)

riverdistmat <- distancematrix.cluster %>%
  dplyr::select(c("hv001.x", "hv001.y", "riverdist")) %>%
  dplyr::rename(distance = riverdist) %>%
  dplyr::mutate(distance = distance/1e3)

#..............................................................
# make distance data
#..............................................................
# ibS
ibS.gcdistmat <- dplyr::left_join(ibS, gcdistmat, by = c("hv001.x", "hv001.y")) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "relatedness", "K1", "K2", "distance")) %>%
  dplyr::mutate(distance = ifelse(K1 == K2, 0, distance))
ibS.roaddistmat <- dplyr::left_join(ibS, roaddistmat, by = c("hv001.x", "hv001.y")) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "relatedness", "K1", "K2", "distance")) %>%
  dplyr::mutate(distance = ifelse(K1 == K2, 0, distance))
ibS.riverdistmat <- dplyr::left_join(ibS, riverdistmat, by = c("hv001.x", "hv001.y")) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "relatedness", "K1", "K2", "distance")) %>%
  dplyr::mutate(distance = ifelse(K1 == K2, 0, distance))

# ibD
ibD.gcdistmat <- dplyr::left_join(ibD, gcdistmat, by = c("hv001.x", "hv001.y")) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "relatedness", "K1", "K2", "distance")) %>%
  dplyr::mutate(distance = ifelse(K1 == K2, 0, distance))
ibD.roaddistmat <- dplyr::left_join(ibD, roaddistmat, by = c("hv001.x", "hv001.y")) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "relatedness", "K1", "K2", "distance")) %>%
  dplyr::mutate(distance = ifelse(K1 == K2, 0, distance))
ibD.riverdistmat <- dplyr::left_join(ibD, riverdistmat, by = c("hv001.x", "hv001.y")) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "relatedness", "K1", "K2", "distance")) %>%
  dplyr::mutate(distance = ifelse(K1 == K2, 0, distance))

#...................
# get scalars
#...................
distscale <- lapply(list(gcdistmat, roaddistmat, riverdistmat),
                    function(x){
                      return(mean(x$distance))
                    })
distscale <- c(distscale, distscale)


#...................
#  Model framework
#...................

modLL <- tibble::tibble(
  name = c("gcdist-ibs", "roaddist-ibs", "riverdist-ibs", "gcdist-ibd", "roaddist-ibd", "riverdist-ibd"),
  distdata = list(ibS.gcdistmat,
                  ibS.roaddistmat,
                  ibS.riverdistmat,
                  ibD.gcdistmat,
                  ibD.roaddistmat,
                  ibD.riverdistmat),
  scalar = distscale
)

#..............................................................
# Get imp
#..............................................................
clsts <- sort( unique(c(ibS.gcdistmat$K1, ibS.gcdistmat$K2)) )

modLL$global_fii <- purrr::map(modLL$distdata, function(dat){
  # nested sapply loop
  f_ii <- sapply(clsts, function(x){
    ret <- mean( dat$relatedness[c(dat$K1 %in% x & dat$K2 %in% x)] )
    return(ret)
  })
  global_fii <- mean(f_ii[!is.nan(f_ii)])
  return(global_fii)
})
# clusters that we should impute on
imputclst <- mtdt %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 1)
imputclst <- imputclst$hv001
# put in model
modLL$imputclst <- lapply(1:nrow(modLL), function(x) return(imputclst))

modLL$ret <- purrr::pmap(modLL, get_distance_geno_likelihood)

modLL$LogLik <- purrr::map(modLL$ret, "LL")
modLL <- modLL %>%
  tidyr::unnest(cols = "LogLik")

#..............................................................
# save out
#..............................................................
modLL.out <- modLL %>%
  dplyr::select(-c("distdata"))
dir.create("results/distance_likelihoods_local/", recursive = T)
saveRDS(modLL.out, file = "results/distance_likelihoods_local/model_likelihoods.RDS")

