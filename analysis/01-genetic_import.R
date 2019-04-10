#----------------------------------------------------------------------------------------------------
# Purpose of this script is simply to filter the genetic data that was previously
# filtered and included in the Big Barcode paper by Bob Verity
#
#----------------------------------------------------------------------------------------------------


# libraries and imports
gdrive <- tcltk::tk_choose.dir()
library(tidyverse)


# Bob Verity filtered this for the
# big barcode paper
mipbigpanel <- readRDS(paste0(gdrive, "/data/raw_data/mip_genetic_data_from_Bob_big_barcode/biallelic_processed.rds"))
drcsmpls <- mipbigpanel$samples$Country == "DRC"

# Subset to DRC Samples
drcmips <- MIPanalyzer::filter_samples(x = mipbigpanel,
                                       sample_filter = drcsmpls,
                                       description = "Subset to DRC Samples")


saveRDS(drcmips, file = paste0(gdrive, "/data/derived_data/cd2013_dhs_bigbarcodemips.rds"))



