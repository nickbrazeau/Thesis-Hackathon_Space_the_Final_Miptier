####################################################################################
## Purpose:
##
## Notes:
####################################################################################
library(tidyverse)
library(raster)
library(tidygraph)
library(ggraph)
source("Space_the_Final_Miptier/R/themes.R")
set.seed(48)
#......................
# read data
#......................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
mtdt <- readRDS("data/derived_data/sample_metadata.rds")
#..............................................................
# Idenity samples with IBD that are greater than half the
# genome (i.e. at or above meiotic siblings)
# These represent likely recent transmission events
#..............................................................
ibD.meiotic <- ibD %>%
  dplyr::filter(malecotf >= 0.5)
ibdmeioticsmpls <- as.character(ibD.meiotic$smpl1)
ibdmeioticsmpls <- c(ibdmeioticsmpls, as.character(ibD.meiotic$smpl2))

#..............................................................
# Coordinate for meiotic pair edges
# need to jitter very slightly to account for within cluster
# edges -- otherwise lines edges can't be drawn
#..............................................................
# coord pts
coords.pts <- mtdt %>%
  dplyr::select(c("name", "hv001", "longnum", "latnum")) %>%
  dplyr::filter(name %in% ibdmeioticsmpls)
# add spatial jitter
coords.pts$long_jitter <- coords.pts$longnum + runif(n = nrow(coords.pts), min = -1e-3, max = 1e-3)
coords.pts$lat_jitter <- coords.pts$latnum + runif(n = nrow(coords.pts), min = -1e-3, max = 1e-3)
coords.pts <- coords.pts %>%
  dplyr::select(c("name", "hv001", "long_jitter", "lat_jitter"))
# long connections
colnames(coords.pts)[1] <- "smpl1"
ibD.meiotic <- dplyr::left_join(ibD.meiotic, coords.pts, by = "smpl1")
colnames(coords.pts)[1] <- "smpl2"
ibD.meiotic <- dplyr::left_join(ibD.meiotic, coords.pts, by = "smpl2")
colnames(coords.pts)[1] <- "name"

# write this out for later use
saveRDS(ibD.meiotic, file = "data/derived_data/bigbarcode_genetic_data/meiotic_sibs_mipanalyzer.DRCibD.long.mtdt.rds")
saveRDS(coords.pts, file = "data/derived_data/meiotic_sibs_coordpts.rds")



#............................................................
# permutations?
#...........................................................

