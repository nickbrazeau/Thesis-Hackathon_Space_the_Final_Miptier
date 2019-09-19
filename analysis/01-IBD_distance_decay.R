#-----------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to look
# at decay of IBD with distance
#-----------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

#....................................................................................
# Import Genetic Data and Subset to the DRC Samples
#....................................................................................
ibd <- readRDS("data/derived_data/bigbarcode_genetic_data/DRCIBD_polarized_biallelic_processed.long.rds")
ibs <- readRDS("data/derived_data/bigbarcode_genetic_data/DRCIBS_polarized_biallelic_processed.long.rds")


#....................................................................................
# Import Distance matrix
#....................................................................................
distancematrix.cluster <- readRDS("data/distance_data/distancematrix_bycluster.rds")
distancematrix.smpl <- readRDS("data/distance_data/distancematrix_bysmpl.rds") %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "gcdistance", "roaddistance"))


#.....................
# EDA of Our Genetic Data
#.....................
# Pairwise Comparisons


#....................................................................................
# IBS Distance Decay -- Smpl Level
#....................................................................................

IBSdist <- dplyr::inner_join(ibs, distancematrix.smpl)




