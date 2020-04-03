#................................................................................................................................
# Purpose of this section is to combine the genetic and mtdt data
# at the CLUSTER level
#................................................................................................................................
library(tidyverse)
source("R/pairwise_helpers.R")
#..............................................................
# Read in genetic dat
#..............................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.rds")
ibS <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS.long.rds")

#..............................................................
# Read in Metadata
#..............................................................
smplmtdt <- readRDS(file = "data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("name", "barcode", "hv001", "adm1name", "latnum", "longnum"))

#..............................................................
# Bring together ibD
#..............................................................
colnames(smplmtdt)[1] <- "smpl1"
ibD <- dplyr::left_join(ibD, smplmtdt, by = "smpl1")
colnames(smplmtdt)[1] <- "smpl2"
ibD <- dplyr::left_join(ibD, smplmtdt, by = "smpl2")

#..............................................................
# Bring together ibS
#..............................................................
colnames(smplmtdt)[1] <- "smpl1"
ibS <- dplyr::left_join(ibS, smplmtdt, by = "smpl1")
colnames(smplmtdt)[1] <- "smpl2"
ibS <- dplyr::left_join(ibS, smplmtdt, by = "smpl2")

#..............................................................
# Save out
#..............................................................
saveRDS(ibD,
        "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")
saveRDS(ibS,
        "data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibS.long.mtdt.rds")




