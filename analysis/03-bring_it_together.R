#----------------------------------------------------------------------------------------------------
# Purpose of this script is simply to collate all the various data sources that we have
# collected
#
#----------------------------------------------------------------------------------------------------

# libraries and imports
gdrive <- tcltk::tk_choose.dir()
library(tidyverse)

#...........................................................
# Epi Data
#...........................................................
dt <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_kids_dhs_recode.rds"))

# merge in Epi
drcmips$samples <- drcmips$samples %>%
  dplyr::rename(sh312 = Barcode) %>%
  dplyr::mutate(sh312 = tolower(sh312)) %>%
  dplyr::left_join(x=., y = dt, by = "sh312")



#...........................................................
# Genetic Data
#...........................................................
# Bob Verity filtered this for the
# big barcode paper
mipbigpanel <- readRDS(paste0(gdrive, "/data/raw_data/mip_genetic_data_from_Bob_big_barcode/biallelic_processed.rds"))

drcsmpls <- mipbigpanel$samples$Country == "DRC"

# Subset to DRC Samples
drcmips <- MIPanalyzer::filter_samples(x = mipbigpanel,
                                       sample_filter = drcsmpls,
                                       description = "Subset to DRC Samples")


# using Bob Verity's MIPAnalyzer package to calculate genetic distance
drcmips$gendistances <- list()
drcmips$gendistances$IBS_verity_maskhets  <- MIPanalyzer::get_IBS_distance(     x = drcmips,
                                                                                ignore_het = T)
drcmips$gendistances$IBS_verity_majallele <- MIPanalyzer::get_IBS_distance(     x = drcmips,
                                                                                ignore_het = T)
drcmips$gendistances$DAB_malariagen_wsaf  <- MIPanalyzer::get_genomic_distance( x = drcmips)

#...........................................................
# Spatial Data
#...........................................................
drc_space <- list()
clst <- dt[,c("hv001", "longnum", "latnum", "travel_times_2015", "adm1fips", "adm1fipsna", "adm1salbna",
                      "adm1salbco", "adm1dhs", "adm1name", "dhsregco", "dhsregna",
                      "urban_rura", "geometry")]

clst <- clst[!duplicated(clst), ]

drc_space$mtdt <- clst
drc_space$spatialdist <- geosphere::distGeo(drc_space$mtdt[,c("longnum", "latnum")])













