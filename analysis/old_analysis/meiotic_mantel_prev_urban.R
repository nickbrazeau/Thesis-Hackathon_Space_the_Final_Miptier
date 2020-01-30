source("R/pairwise_helpers.R")
library(tidyverse)
mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("name", "barcode", "hv001", "longnum", "latnum"))

ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::select(c("hv001", "urban_rura")) %>%
  dplyr::mutate(urban_rura_ext = ifelse(urban_rura == "R", 10000, 2000))

mtdt <- dplyr::left_join(mtdt, ge, by = "hv001")


ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds") %>%
  dplyr::filter(malecotf >= 0.5) %>%
  dplyr::filter(hv001.x != hv001.y)

#..............................................................
# get urban
#..............................................................

urban <- raster::raster("data/derived_data/urbanicity_raster/urbanicity.grd")
urbanmean <- rep(NA, nrow(mtdt))
for (i in 1:nrow(mtdt)){
  urbanmean[i] <- raster::extract(x = urban,
                                  y = sf::as_Spatial(mtdt$geometry[i]),
                                  buffer = mtdt$urban_rura_ext[i],
                                  fun = mean,
                                  na.rm = T)
}
mtdt$urban <- urbanmean

#..............................................................
# diff
#..............................................................
urban.pts <- mtdt %>%
  dplyr::select(c("name", "urban")) %>%
  dplyr::rename(smpl1 = name,
                urbanmean = urban)

ibD.urban <- dplyr::left_join(ibD, urban.pts, by = "smpl1")


colnames(urban.pts)[1] <- "smpl2"
ibD.urban <- dplyr::left_join(ibD.urban, urban.pts, by = "smpl2")

ibD.urban <- ibD.urban %>%
  dplyr::mutate(
    urban.euclidean = sqrt( (urbanmean.x - urbanmean.y)^2 )
  )

#..............................................................
# Prevalence
#..............................................................

prev <- raster::raster("data/derived_data/MAPrasters/parasiterate.grd")
prevmean <- rep(NA, nrow(mtdt))
mtdt$urban_rura_ext[mtdt$hv001 == 54] <- 6000 # missing data otherwise
for (i in 1:nrow(mtdt)){
  prevmean[i] <- raster::extract(x = prev,
                                 y = sf::as_Spatial(mtdt$geometry[i]),
                                 buffer = mtdt$urban_rura_ext[i],
                                 fun = mean,
                                 na.rm = T)
}
mtdt$prev <- prevmean

#..............................................................
# diff
#..............................................................
prev.pts <- mtdt %>%
  dplyr::select(c("name", "prev")) %>%
  dplyr::rename(smpl1 = name,
                prevmean = prev)

ibD.urban.prev <- dplyr::left_join(ibD.urban, prev.pts, by = "smpl1")


colnames(prev.pts)[1] <- "smpl2"
ibD.urban.prev <- dplyr::left_join(ibD.urban.prev, prev.pts, by = "smpl2")

ibD.urban.prev <- ibD.urban.prev %>%
  dplyr::mutate(
    prev.euclidean = sqrt( (prevmean.x - prevmean.y)^2 )
  )

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

# covar distance
urban.wide <- make_symmat(ibD.urban.prev, urban.euclidean)
prev.wide <- make_symmat(ibD.urban.prev, prev.euclidean)

# genetic distance
ibD.wide <- make_symmat(ibD.urban.prev, "malecotf")



#...........................................................
# Perform Mantel Tests
#...........................................................
mantel_testsdf <- tibble::tibble(comparison = c("urban", "prev"),
                                 covardist = list(urban.wide, prev.wide))
mantel_testsdf$ibD <- lapply(1:nrow(mantel_testsdf), function(x) return(ibD.wide))



mantel_testsdf$manteltest_ibD <- purrr::pmap(mantel_testsdf[,c("ibD", "covardist")],
                                             function(ibD, covardist){
                                               ret <- ape::mantel.test(m1 = ibD,
                                                                       m2 = covardist,
                                                                       nperm = 1e4,
                                                                       graph = T,
                                                                       alternative = "two.sided")
                                             })

out <- mantel_testsdf %>%
  dplyr::select(c("comparison", dplyr::starts_with("mantel")))
dir.create(path = "~/Documents/GitHub/Space_the_Final_Miptier/results/simresults/", recursive = T)
saveRDS(out, "~/Documents/GitHub/Space_the_Final_Miptier/results/simresults/mantel_testsdf_covars_meiotic.RDS")







