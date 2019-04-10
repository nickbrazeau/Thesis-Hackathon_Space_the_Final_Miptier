#----------------------------------------------------------------------------------------------------
# Purpose of this script is simply to collate all the various data sources that we have
# collected
#
#----------------------------------------------------------------------------------------------------

# libraries and imports
gdrive <- tcltk::tk_choose.dir()
library(tidyverse)



#...........................................................
# Read in Genetic Data
#...........................................................
drcmips <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_dhs_bigbarcodemips.rds"))


#...........................................................
# Read in Epi Data
#...........................................................
dt <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_kids_dhs_recode.rds"))

# merge in Epi
drcmips$samples <- drcmips$samples %>%
  dplyr::rename(sh312 = Barcode) %>%
  dplyr::mutate(sh312 = tolower(sh312)) %>%
  dplyr::left_join(x=., y = dt, by = "sh312")




#...........................................................
# Calculate Genetic Distances
#...........................................................
# using Bob Verity's MIPAnalyzer package to calculate genetic distance
drcmips$gendistances <- list()
drcmips$gendistances$IBS_verity_maskhets  <- MIPanalyzer::get_IBS_distance(     x = drcmips,
                                                                                ignore_het = T)
drcmips$gendistances$IBS_verity_majallele <- MIPanalyzer::get_IBS_distance(     x = drcmips,
                                                                                ignore_het = F)
drcmips$gendistances$DAB_malariagen_wsaf  <- MIPanalyzer::get_genomic_distance( x = drcmips)





#....................
# Aggregate by admin level
#.....................
type = c("mean", "median", "trimmed")

smpl_admin_key1 <- drcmips$samples %>%
  dplyr::select(c("sh312", "ADM1NAME")) %>%
  dplyr::mutate( item1 = as.character( seq(1, nrow(.)) )) %>% # key
  dplyr::rename(smpl1 = sh312,
                admin1_1 = ADM1NAME)


smpl_admin_key2 <- drcmips$samples %>%
  dplyr::select(c("sh312", "ADM1NAME")) %>%
  dplyr::mutate( item2 = as.character( seq(1, nrow(.)) )) %>% # key
  dplyr::rename(smpl2 = sh312,
                admin1_2 = ADM1NAME)


IBS_verity_maskhets.tidy <-
  broom::tidy(as.dist( t( drcmips$gendistances$IBS_verity_maskhets) ))

IBS_verity_maskhets.tidy <- IBS_verity_maskhets.tidy %>%
  dplyr::left_join(x=., y =smpl_admin_key1, by = "item1")

IBS_verity_maskhets.tidy <- IBS_verity_maskhets.tidy %>%
  dplyr::left_join(x=., y =smpl_admin_key2, by = "item2")


ret <- IBS_verity_maskhets.tidy %>%
  dplyr::group_by(admin1_1, admin1_2) %>%
  dplyr::summarise(
    adim_dist = mean(distance, na.rm = T)
  )




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













