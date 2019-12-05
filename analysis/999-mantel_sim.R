#..............................................................
# Mantel Simulations for spatial and genetic auotocorrelation
#..............................................................
library(tidyverse)
source("~/Documents/GitHub/Space_the_Final_Miptier/R/pairwise_helpers.R")
set.seed(48)
#..............................................................
# get metadata
#..............................................................
drcsmpls <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/drcsmpls_foruse.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("id", "hv001"))

#....................................................................................
# Import Genetic Data
#....................................................................................
ibD <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD_polarized_biallelic_processed.long.rds")
ibS <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS_polarized_biallelic_processed.long.rds")

#.....................
# Manipulate for Decay Comparisons
#.....................
# merge in cluster informations
colnames(drcsmpls)[1] <- "smpl1"
ibD <- dplyr::left_join(ibD, drcsmpls, by = "smpl1")
ibS <- dplyr::left_join(ibS, drcsmpls, by = "smpl1")

colnames(drcsmpls)[1] <- "smpl2"
ibD <- dplyr::left_join(ibD, drcsmpls, by = "smpl2")
ibS <- dplyr::left_join(ibS, drcsmpls, by = "smpl2")

#.......................................................................
# Import Distance matrix
#.......................................................................
distancematrix.cluster <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/distancematrix_bycluster.rds")

# IBS/D Mantel Tests
# The mantel test assess the correlation between a genetic distance and a geographical distance. In this case, we will use IBS and IBD as the measures of genetic distance (although $F_{st}$, Nei's $\pi$, etc. could be used). Similarly, geographical distance is measured with either greater circle, road, or river distance. Another thank you to  Emmanuel Pardis for making this analyses easy with ape. We will use 10,000 permutations. _Note, we will use two-sided p-values_.

# need to make the distance matrices "wide"
ibSdist <- long_distance_matrix_join(x=ibS, y=distancematrix.cluster,
                                     by = c("hv001.x", "hv001.y")) %>%
  dplyr::mutate(gcdistance = gcdistance/1e3,
                roaddistance = roaddistance/1e3,
                riverdist = riverdist/1e3,
                riverdist = as.numeric(riverdist)) # remove units
# take care of cluster diagnonals
ibSdist$gcdistance[ibSdist$hv001.x == ibSdist$hv001.y] <- 0
ibSdist$roaddistance[ibSdist$hv001.x == ibSdist$hv001.y] <- 0
ibSdist$riverdist[ibSdist$hv001.x == ibSdist$hv001.y] <- 0


ibDdist <- long_distance_matrix_join(x=ibD, y=distancematrix.cluster,
                                     by = c("hv001.x", "hv001.y")) %>%
  dplyr::mutate(gcdistance = gcdistance/1e3,
                roaddistance = roaddistance/1e3,
                riverdist = riverdist/1e3,
                riverdist = as.numeric(riverdist)) # remove units


# take care of cluster diagnonals
ibDdist$gcdistance[ibDdist$hv001.x == ibDdist$hv001.y] <- 0
ibDdist$roaddistance[ibDdist$hv001.x == ibDdist$hv001.y] <- 0
ibDdist$riverdist[ibDdist$hv001.x == ibDdist$hv001.y] <- 0


#...........................................................
# Set up Matrices Mantel Tests
#...........................................................
# take these to wide format for distance matrices
make_symmat <- function(dist, var){
  var <- enquo(var)
  dist.wide <- dist %>%
    dplyr::select(c("smpl1", "smpl2", !!var)) %>%
    tidyr::spread(., key = "smpl1", value = !!var)
  # need to make this matrix full and have a diagonal
  dist.wide[,1] <- NA # these were sample names
  dist.wide <- as.matrix( rbind(dist.wide, NA) )
  diag(dist.wide) <- 0
  # make symmetric matrix
  dist.wide[ lower.tri(dist.wide) ] <- dist.wide[ upper.tri(dist.wide) ]

  return(dist.wide)
}

# geographic distance
gc.wide <- make_symmat(ibSdist, gcdistance)
road.wide <- make_symmat(ibSdist, roaddistance)
river.wide <- make_symmat(ibSdist, riverdist)

# genetic distance
ibS.wide <- make_symmat(ibSdist, "hammings")
ibD.wide <- make_symmat(ibDdist, "malecotf")


#...........................................................
# Perform Mantel Tests
#...........................................................
mantel_testsdf <- tibble::tibble(comparison = c("gc", "road", "river"),
                                 geographicdist = list(gc.wide, road.wide, river.wide))
mantel_testsdf$ibS <- lapply(1:nrow(mantel_testsdf), function(x) return(ibS.wide))
mantel_testsdf$ibD <- lapply(1:nrow(mantel_testsdf), function(x) return(ibD.wide))

mantel_testsdf$manteltest_ibS <- purrr::pmap(mantel_testsdf[,c("ibS", "geographicdist")],
                                             function(ibS, geographicdist){
                                               ret <- ape::mantel.test(m1 = ibS,
                                                                       m2 = geographicdist,
                                                                       nperm = 1e5,
                                                                       graph = T,
                                                                       alternative = "two.sided")
                                             })


mantel_testsdf$manteltest_ibD <- purrr::pmap(mantel_testsdf[,c("ibD", "geographicdist")],
                                             function(ibD, geographicdist){
                                               ret <- ape::mantel.test(m1 = ibD,
                                                                       m2 = geographicdist,
                                                                       nperm = 1e4,
                                                                       graph = T,
                                                                       alternative = "two.sided")
                                             })

out <- mantel_testsdf %>%
  dplyr::select(c("comparison", dplyr::starts_with("mantel")))
dir.create(path = "~/Documents/GitHub/Space_the_Final_Miptier/results/simresults/", recursive = T)
saveRDS(out, "~/Documents/GitHub/Space_the_Final_Miptier/results/simresults/mantel_testsdf.RDS")


