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
ibD <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
# rename for recursive function
ibD <- ibD %>%
  dplyr::rename(relatedness = malecotf) %>%
  dplyr::select(c("smpl1", "smpl2", "relatedness", "hv001.x", "hv001.y"))

# ibs
ibS <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS.long.mtdt.rds")
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
  fulldata = list(ibS.gcdistmat,
                  ibS.roaddistmat,
                  ibS.riverdistmat,
                  ibD.gcdistmat,
                  ibD.roaddistmat,
                  ibD.riverdistmat),
  scalar = distscale
)


#..............................................................
# Massive parallelization
#..............................................................
clsts <- sort( unique(c(ibS.gcdistmat$K1, ibS.gcdistmat$K2)) )
clsts.combns <- as.data.frame( t( combn(clsts, 2) ) )
colnames(clsts.combns) <- c("K1sub", "K2sub")


modLL$distdata <- purrr::map(modLL$fulldata, function(dat){
  # nested lapply to subset to specific row of data that we want
  out <- apply(clsts.combns, 1, function(combns){
    ret <- dat %>%
      dplyr::filter(c(K1 %in% as.vector(combns) | K2 %in% as.vector(combns)))
    return(ret)
  })
  return(out)
})

# subset and unnest
modLL$Kselect <- lapply(1:nrow(modLL), function(x) return(clsts.combns))
modLL <- modLL %>%
  dplyr::select(-c("fulldata")) %>%
  tidyr::unnest(cols = c(Kselect, distdata))

#..............................................................
# Save out here so we don't have to read distdata all back in
# and create mem bottleneck for rmd downstream
#..............................................................
easyout <- modLL %>%
  dplyr::select(-c("distdata"))
saveRDS(object = easyout, file = "data/derived_data/Likelihood_backend_pairwise_comparisons_params_nodata.RDS")


#..............................................................
# Add in data for imputation on clusters with only one
# sample
#..............................................................
# global within cluster ibd to impute for clusters with only one sample
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

# for slurm on LL
dir.create("results/distance_likelihoods", recursive = T)
setwd("results/distance_likelihoods/")
ntry <- 1028 # max number of nodes allowed
sjob <- rslurm::slurm_apply(f = get_distance_geno_likelihood,
                            params = modLL,
                            jobname = 'distance_likelihoods_recursion',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = 256000,
                                                 array = sprintf("0-%d%%%d",
                                                                 ntry,
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "10:00:00"))



