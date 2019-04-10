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
source("R/00-gendist_liftover_functions.R")
drcmips$results <- list()
drcmips$results$admin_gen_dist <- list()
# Verity IBS mask all hets
drcmips$results$admin_gen_dist$IBS_verity_maskhets <- getadmin_gendist_summary(
                                                                  mipanalyzerobject_samples = drcmips$samples,
                                                                  mipanalyzerobject_gendistmat = drcmips$gendistances$IBS_verity_maskhets,
                                                                  type = "mean" # TODO fix matcharg
                                                                  )

# Verity IBS force all hets
drcmips$results$admin_gen_dist$IBS_verity_majallele <- getadmin_gendist_summary(
  mipanalyzerobject_samples = drcmips$samples,
  mipanalyzerobject_gendistmat = drcmips$gendistances$IBS_verity_majallele,
  type = "mean" # TODO fix matcharg
)

# MalariaGEN DAB
drcmips$results$admin_gen_dist$DAB_malariagen_wsaf <- getadmin_gendist_summary(
  mipanalyzerobject_samples = drcmips$samples,
  mipanalyzerobject_gendistmat = drcmips$gendistances$DAB_malariagen_wsaf,
  type = "mean" # TODO fix matcharg
)






#...........................................................
# Spatial Data
#...........................................................
drcmips$spatialdist <- list()

clst <- dt[,c("hv001", "longnum", "latnum", "travel_times_2015", "adm1fips", "adm1fipsna", "adm1salbna",
                      "adm1salbco", "adm1dhs", "adm1name", "dhsregco", "dhsregna",
                      "urban_rura", "geometry")]
clst <- clst[!duplicated(clst), ]

admin_centroids_clst_wght <- clst %>%
  dplyr::group_by(adm1name) %>%
  dplyr::summarise(
    admin_long = mean(longnum),
    admin_lat = mean(latnum)
  )
sf::st_geometry(admin_centroids_clst_wght) <- NULL


# Greater Circle Distance
drcmips$spatialdist$admin_gc_dist <- MIPanalyzer::get_spatial_distance(lat = unlist( admin_centroids_clst_wght[,"admin_lat"] ),
                                  long = unlist( admin_centroids_clst_wght[,"admin_long"]) )

## make key for liftover
adminkey_admin1 <- data.frame(
  item1 = as.character( seq(1, 26) ),
  admin1_1 = levels(factor(admin_centroids_clst_wght$adm1name))
)

adminkey_admin2 <- data.frame(
  item2 = as.character( seq(1, 26) ),
  admin1_2 = levels(factor(admin_centroids_clst_wght$adm1name))
)

# Save spatial distance out as a result
admin_gc_dist.tidy <- broom::tidy(as.dist(drcmips$spatialdist$admin_gc_dist)) %>%
  dplyr::mutate(item1 = gsub("admin_long", "", item1),
                item2 = gsub("admin_long", "", item2))

admin_gc_dist.tidy <- dplyr::left_join(admin_gc_dist.tidy, adminkey_admin1, by = "item1")
admin_gc_dist.tidy <- dplyr::left_join(admin_gc_dist.tidy, adminkey_admin2, by = "item2")

# save out
drcmips$results$admin_gc_dist <- admin_gc_dist.tidy





saveRDS(drcmips, paste0(gdrive, "/data/derived_data/cd2013_gen_space_epi_final.rds"))

