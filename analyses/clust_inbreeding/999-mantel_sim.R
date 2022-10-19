#..............................................................
# Mantel Simulations for spatial and genetic auotocorrelation
#..............................................................
library(tidyverse)
source("R/pairwise_helpers.R")
set.seed(48)

#....................................................................................
# Import Genetic Data
#....................................................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds") %>%
  dplyr::mutate(hv001.x = as.character(hv001.x),
                hv001.y = as.character(hv001.y))

#.......................................................................
# Import Distance matrix
#.......................................................................
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")

# IBD Mantel Tests
# The mantel test assess the correlation between a genetic distance and a geographical distance. In this case, we will use IBS and IBD as the measures of genetic distance (although $F_{st}$, Nei's $\pi$, etc. could be used). Similarly, geographical distance is measured with either greater circle, road, or river distance. Another thank you to  Emmanuel Pardis for making this analyses easy with ape. We will use 10,000 permutations. _Note, we will use two-sided p-values_z
ibDdist <- long_distance_matrix_join(x=ibD, y=distancematrix.cluster,
                                     by = c("hv001.x", "hv001.y")) %>%
  dplyr::mutate(gcdistance = gcdistance/1e3,
                roaddistance = roaddistance/1e3,
                airportdistance = airportdistance/1e3)

# take care of cluster diagnonals
ibDdist$gcdistance[ibDdist$hv001.x == ibDdist$hv001.y] <- 0
ibDdist$roaddistance[ibDdist$hv001.x == ibDdist$hv001.y] <- 0
ibDdist$airportdistance[ibDdist$hv001.x == ibDdist$hv001.y] <- 0


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
  dist.wide[ lower.tri(dist.wide) ] <- t(dist.wide)[ upper.tri(dist.wide) ]

  return(dist.wide)
}

# geographic distance
gc.wide <- make_symmat(ibDdist, gcdistance)
road.wide <- make_symmat(ibDdist, roaddistance)
airport.wide <- make_symmat(ibDdist, airportdistance)

# genetic distance
ibD.wide <- make_symmat(ibDdist, "malecotf")


#...........................................................
# Perform Mantel Tests
#...........................................................
mantel_testsdf <- tibble::tibble(comparison = c("gc", "road", "airport"),
                                 geographicdist = list(gc.wide, road.wide, airport.wide))
mantel_testsdf$ibD <- lapply(1:nrow(mantel_testsdf), function(x) return(ibD.wide))

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
dir.create(path = "results/mantel_perm_results/", recursive = T)
saveRDS(out, "results/mantel_perm_results/mantel_testsdf.RDS")


