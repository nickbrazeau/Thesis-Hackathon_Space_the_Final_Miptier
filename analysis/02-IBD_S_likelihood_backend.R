#..............................................................
# Purpose of this script is to calculate distance
# likelihood
#..............................................................
library(tidyverse)
source("R/pairwise_helpers.R")
#..............................................................
# Function
#..............................................................

get_distance_geno_likelihood <- function(name, relatedness, distmat, clsts){
  #..............................................................
  # functions
  #..............................................................
  gaussdist <- function(x){
    x <- x/1e3 # m to km
    ret <- (1 / sqrt(2*pi) ) * exp(-x)
    return(ret)
  }

  #..............................................................
  # setup
  #..............................................................
  clsts <- unique(clsts)
  fuu <- sapply(clsts, function(x){
    ret <- mean( relatedness$relatedness[c(relatedness$hv001.x %in% x | relatedness$hv001.y %in% x)] )
    return(ret)
  })


  for (i in 1:length(clsts)) {
    cat("Iteration i: ", i)
    for (j in 1:length(clsts)) {

      # init
      LL <- 0

      # first term
      frsttrm <- 0
      for (u in 1:length(clsts)) {
        if ( clsts[u] == clsts[i] ){
          diu <- 0
        } else {
          diu <- gaussdist( distmat$distance[distmat$hv001.x == clsts[i] & distmat$hv001.y == clsts[u]] )
        }
        if ( clsts[u] == clsts[j] ){
          dju <- 0
        } else {
          dju <- gaussdist( distmat$distance[distmat$hv001.x == clsts[j] & distmat$hv001.y == clsts[u]] )
        }

        frsttrm <- frsttrm + (fuu[u] * diu * dju)

      }

      # second term
      scndterm <- 0
      # now loop through
      for (u in 1:length(clsts)) {
        for (v in length(clsts):1) { # go backwards so not same deme
          fuv <- mean( relatedness$relatedness[ relatedness$hv001.x %in% c( clsts[u], clsts[v] ) ] ) # mean within demes relatedness
          if (clsts[u] == clsts[i]) {
            diu <- 0
          } else {
            diu <- gaussdist( distmat$distance[distmat$hv001.x == clsts[i] & distmat$hv001.y == clsts[u]] )
          }

          if (clsts[v] == clsts[j]) {
            djv <- 0
          } else {
            djv <- gaussdist( distmat$distance[distmat$hv001.x == clsts[j] & distmat$hv001.y == clsts[v]] )
          }
          scndterm <- scndterm + (fuv*diu*djv)
        }
      }

      # add to LL
      LL <- LL + sum(log(frsttrm), log(scndterm))
    } # end loop j
  } # end loop i
  return(LL)
}


#..............................................................
# Model Setup
#..............................................................
#...................
# meta Data
#...................
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

# ibs
ibS <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS_polarized_biallelic_processed.long.rds")
# merge in cluster informations
colnames(drcsmpls)[1] <- "smpl1"
ibS <- dplyr::left_join(ibS, drcsmpls, by = "smpl1")

colnames(drcsmpls)[1] <- "smpl2"
ibS <- dplyr::left_join(ibS, drcsmpls, by = "smpl2")
# rename for recursive function
ibS <- ibS %>%
  dplyr::rename(relatedness = hammings)

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

roaddistmat <- distancematrix.cluster %>%
  dplyr::select(c("hv001.x", "hv001.y", "roaddistance")) %>%
  dplyr::rename(distance = roaddistance) %>%
  dplyr::mutate(distance = distance/1e3)


riverdistmat <- distancematrix.cluster %>%
  dplyr::select(c("hv001.x", "hv001.y", "riverdist")) %>%
  dplyr::rename(distance = riverdist) %>%
  dplyr::mutate(distance = distance/1e3)



#...................
#  Model framework
#...................

modLL <- tibble::tibble(
  name = c("gcdist-ibs", "roaddist-ibs", "riverdist-ibs", "gcdist-ibd", "roaddistibd", "riverdist-ibd"),
  relatedness = list(ibS, ibS, ibS, ibD, ibD, ibD),
  distmat = list(gcdistmat, roaddistmat, riverdistmat, gcdistmat, roaddistmat, riverdistmat),
  clsts = list(drcsmpls$hv001)
)


nick = purrr::pmap(modLL, get_distance_geno_likelihood)


# for slurm on LL
dir.create("results/distance_likelihoods", recursive = T)
setwd("results/distance_likelihoods/")
ntry <- nrow(modLL)
sjob <- rslurm::slurm_apply(f = get_distance_geno_likelihood,
                            params = modLL,
                            jobname = 'distance_likelihoods_recursion',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = 64000,
                                                 array = sprintf("0-%d%%%d",
                                                                 ntry,
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))



