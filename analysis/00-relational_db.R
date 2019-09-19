#-----------------------------------------------------------------------------------------------------------------------------------------
# Purpose of this script is to set up relational data sets
#-----------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

drcgeo <- readRDS("data/derived_data/DRCgeo_points.rds")
covarmatrix <- readRDS("data/derived_data/covariatematrix.rds")
distancematrix <- readRDS("data/derived_data/distancematrix.rds")

# genetic data
mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.)))
coi <- readRDS("data/derived_data/coi_data.rds")

ibd <- readr::read_tsv("data/derived_data/bigbarcode_genetic_data/IBD_polarized_biallelic_processed.long.tab.txt", col_names = F) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "nsites", "malecotf"))

ibs <- read.table("data/derived_data/bigbarcode_genetic_data/IBS_polarized_biallelic_processed.long.tab.txt",
                  header = T, stringsAsFactors = F)
rownames(ibs) <- unlist(ibs[,1])
ibs <- ibs[,-1]
ibs <- broom::tidy( as.dist(ibs) ) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "hammingibs"))
