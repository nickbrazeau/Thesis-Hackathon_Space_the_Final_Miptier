#..............................................................
# Moran Simulations for spatial auotocorrelation of centrality
#..............................................................
library(tidyverse)
source("~/Documents/GitHub/Space_the_Final_Miptier/R/pairwise_helpers.R")
set.seed(48)
#..............................................................
# Import metadata
#..............................................................
drcsmpls <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/drcsmpls_foruse.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("id", "hv001")) %>%
  dplyr::rename(name = id)

mtdt <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(name = id) %>%
  dplyr::select(c("name", "country", "hv001", "adm1name", "longnum", "latnum")) %>%
  dplyr::filter(name %in% drcsmpls$name)
#.......................................................................
# Import Distance matrix
#.......................................................................
distancematrix.cluster <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/distancematrix_bycluster.rds")

#....................................................................................
# Import Centrality Measures
#....................................................................................
cent <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/centrality_measures/centrality_final_wi_and_knn5.RDS")
# need these for join
cent.wide <- cent %>%
  dplyr::left_join(x = ., y = mtdt, by = "name") %>%
  dplyr::select(c("name", dplyr::starts_with("centrality"), dplyr::starts_with("hv001")))

#.......................................................................
# Import genetic data
#.......................................................................
# need to make the distance matrices "wide"
# Don't need IBS calculates here but need to replicate duplicated coordinates
ibS <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS_polarized_biallelic_processed.long.rds")
#.....................
# Manipulate for Comparisons
#.....................
# merge in cluster informations
colnames(drcsmpls)[1] <- "smpl1"
ibS <- dplyr::left_join(ibS, drcsmpls, by = "smpl1")

colnames(drcsmpls)[1] <- "smpl2"
ibS <- dplyr::left_join(ibS, drcsmpls, by = "smpl2")

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

geodist <- ibSdist %>%
  dplyr::select(-c("hammings", "hv001.x", "hv001.y"))

#...........................................................
# Set up Matrices Moran I Tests
#...........................................................
# take these to wide format for distance matrices
make_symmat <- function(dist, var){
  var <- enquo(var)
  dist.wide <- dist %>%
    dplyr::select(c("smpl1", "smpl2", !!var)) %>%
    tidyr::spread(., key = "smpl1", value = !!var)

  names.rows <- dist.wide$smpl2
  names.col <- colnames(dist.wide)
  mssnglstsmpl <- names.col[ ! names.rows %in% c(names.col) ]
  mssnglstsmpl <- mssnglstsmpl[ mssnglstsmpl != "smpl2"]
  names.rows <- c(names.rows, mssnglstsmpl)
  # need to make this matrix full and have a diagonal
  dist.wide[,1] <- NA # these were sample names
  dist.wide <- as.matrix( rbind(dist.wide, NA) )
  diag(dist.wide) <- 0
  # make symmetric matrix
  dist.wide[ lower.tri(dist.wide) ] <- dist.wide[ upper.tri(dist.wide) ]

  rownames(dist.wide) <- names.rows
  colnames(dist.wide) <- NULL
  return(dist.wide)
}

# geographic distance wide
gc.wide <- make_symmat(geodist, gcdistance)
road.wide <- make_symmat(geodist, roaddistance)
river.wide <- make_symmat(geodist, riverdist)


moran_mc_wrapper <- function(datavect, distmat){
  invdist <- 1/distmat
  diag(invdist) <- 0

  ret <- spdep::moran.mc(x = datavect,
                         listw = spdep::mat2listw(invdist),
                         alternative = "greater",
                         nsim = 1e4,
                         spChk = T)
  return(ret)
}


#........................
# Moran's I Calcs for Weighted
#........................
moransIdf <- tibble::tibble(
  cat = c("gc", "road", "water"),
  distmat = list(gc.wide, road.wide, river.wide)
)

# get cent measures
cent <- cent %>%
  dplyr::arrange(name)
wi_datavect <- cent$centrality_pr_wi
names(wi_datavect) <- cent$name
pr_datavect <- cent$centrality_pr_knn5
names(pr_datavect) <- cent$name

# put in cent measure
moransIdf$datavect <- lapply(1:nrow(moransIdf), function(x) return(wi_datavect))
moransIdf$moranIret <- purrr::pmap(moransIdf[,c("distmat", "datavect")], moran_mc_wrapper)

#..............................................................
# Save out
#..............................................................
out <- moransIdf %>%
  dplyr::select(c("cat", "moranIret"))

dir.create(path = "results/simresults/", recursive = T)
saveRDS(out, "results/simresults/moransI_testsdf.RDS")






