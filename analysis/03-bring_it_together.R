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
mipbigpanel <- readRDS(paste0(gdrive, "/data/derived_data/all_bigbarcodemips.rds"))
drcmips <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_dhs_bigbarcodemips.rds"))
# ghanamips <- readRDS(paste0(gdrive, "/data/derived_data/ghana_bigbarcodemips.rds"))
# ghana_drcmips <- readRDS(paste0(gdrive, "/data/derived_data/ghana_drc_bigbarcodemips.rds"))

#...........................................................
# Read in Epi Data
#...........................................................
dt <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_kids_dhs_recode.rds"))

# individual level covariates with underscores in that change on a per cluster basis
nms <- c();for(i in unique(dt$hv001)){
  nms <- c(nms,names(which(apply(dt[dt$hv001==i,grep("_|hv001",names(dt))],2,function(x) length(unique(x)))>1)))
}
nms <- unique(nms)

#...........................................................
# Merge in Epi Data
#...........................................................

# first merge by raw dhs covariates using hhid and hvidx - this is because the barcodes from the samples
# are found both within sh312 and ha62 and ha63. The Barcode variables in drcmips$samples is created from
# all of this information and then we added in hhid and hvidx so we can always reference them.

# match to everything in this way that does not have an _ except for the nms vars
drcmips$samples <- drcmips$samples %>%
  dplyr::left_join(x=., y = dt[,-grep("_",names(dt))[!grep("_",names(dt)) %in% which(names(dt) %in% nms)]],
                   by = c("hhid","hvidx","hv001"))

# second merge by the cluster for all vars that have a _ in (apart from nms) as these are cluster relavant
# and this way we can bring in covariates for the missing children samples at the cluster level
drcmips$samples <- drcmips$samples %>%
  dplyr::left_join(x=., y = dt[match(unique(dt$hv001),dt$hv001),-which(names(dt) %in% nms)] %>%
                     sf::st_set_geometry(NULL) %>%
                     select(matches("_|hv001")),
                   by = "hv001")

#...........................................................
# Add Spatial Distance Data
#...........................................................
drcmips$spatial_distances <- list()

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
sf::st_geometry(clst) <- NULL


# Greater Circle Distance for admin cluster and individuals
# drcmips$spatial_distances$admin_gc_dist <- MIPanalyzer::get_spatial_distance(lat = unlist( admin_centroids_clst_wght[,"admin_lat"] ),
#                                                                        long = unlist( admin_centroids_clst_wght[,"admin_long"]) )
drcmips$spatial_distances$gc_dist <- MIPanalyzer::get_spatial_distance(lat = drcmips$samples$LATNUM,long = drcmips$samples$LONGNUM)
mipbigpanel$spatial_distances$gc_dist <- MIPanalyzer::get_spatial_distance(lat = mipbigpanel$samples$LATNUM,long = mipbigpanel$samples$LONGNUM)


source("R/00-distance_matrix_munging.R")

#...........................................................
# Add Genetic Distance Data
#...........................................................

drcmips <- add_genetic_distances(drcmips)
mipbigpanel <- add_genetic_distances(mipbigpanel)

#...........................................................
# Add Summary Data Frames
#...........................................................

drcmips <- add_dist_summaries(drcmips, select = c("ADM1NAME","hv001","Barcode"), distances = drcmips$genetic_distances)
# drcmips <- add_dist_summaries(drcmips, select = "hv001",distances = drcmips$genetic_distances)
# drcmips <- add_dist_summaries(drcmips, select = "Barcode",distances = drcmips$genetic_distances)
# drcmips <- add_dist_summaries(drcmips, select = "DHSREGNA",distances = drcmips$genetic_distances)

drcmips <- add_dist_summaries(drcmips, select = c("ADM1NAME","hv001","Barcode"), distances = drcmips$spatial_distances)
# drcmips <- add_dist_summaries(drcmips, select = "hv001",distances = drcmips$spatial_distances)
# drcmips <- add_dist_summaries(drcmips, select = "Barcode",distances = drcmips$spatial_distances)
# drcmips <- add_dist_summaries(drcmips, select = "DHSREGNA",distances = drcmips$spatial_distances)

mipbigpanel <- add_dist_summaries(mipbigpanel, select = "Barcode",distances = mipbigpanel$genetic_distances)
mipbigpanel <- add_dist_summaries(mipbigpanel, select = "Barcode",distances = mipbigpanel$spatial_distances)



# Save spatial distance out as a result
saveRDS(drcmips, paste0(gdrive, "/data/derived_data/cd2013_gen_space_epi_final.rds"))
saveRDS(mipbigpanel, paste0(gdrive, "/data/derived_data/all_bigbarcode_gen_space_epi_final.rds"))



###
## Deprecated


# #...........................................................
# # Calculate Genetic Distances
# #...........................................................
# # using Bob Verity's MIPAnalyzer package to calculate genetic distance
# drcmips$gendistances <- list()
# ghana_drcmips$gendistances <- list()
#
# # all end up as either 1 or NA. Meaningless
# drcmips$gendistances$IBS_verity_maskhets  <- MIPanalyzer::get_IBS_distance(x = drcmips,ignore_het = T)
# ghana_drcmips$gendistances$IBS_verity_maskhets  <- MIPanalyzer::get_IBS_distance(x = ghana_drcmips,ignore_het = T)
#
# # works out fine but forces everything to be similar by forcing hets
# drcmips$gendistances$IBS_verity_majallele <- MIPanalyzer::get_IBS_distance(x = drcmips,ignore_het = F)
# ghana_drcmips$gendistances$IBS_verity_majallele  <- MIPanalyzer::get_IBS_distance(x = ghana_drcmips,ignore_het = F)
#
# # dab distance
# drcmips$gendistances$DAB_malariagen_wsaf  <- MIPanalyzer::get_genomic_distance(x = drcmips)
# ghana_drcmips$gendistances$DAB_malariagen_wsaf  <- MIPanalyzer::get_genomic_distance(x = ghana_drcmips)
#
# # ibd distance
# drcmips$gendistances$MLE_IBD  <- MIPanalyzer::inbreeding_mle(x = drcmips,f = seq(0.01,0.99,0.01), ignore_het = F)
# ghana_drcmips$gendistances$MLE_IBD  <- MIPanalyzer::inbreeding_mle(x = ghana_drcmips,f = seq(0.01,0.99,0.01), ignore_het = F)
#
# # mixture distance
# drcmips$gendistances$IBM_0  <- MIPanalyzer::get_IB_mixture(x = drcmips)
# drcmips$gendistances$IBM_05  <- MIPanalyzer::get_IB_mixture(x = drcmips,tol = 0.05)
# ghana_drcmips$gendistances$IBM_05  <- MIPanalyzer::get_IB_mixture(x = ghana_drcmips,tol = 0.05)


#....................
# Aggregate by admin level
#.....................
# source("R/00-gendist_liftover_functions.R")
# drcmips$results <- list()
# drcmips$results$admin_gen_dist <- list()
#
# # Verity IBS mask all hets
# drcmips$results$admin_gen_dist$IBS_verity_maskhets <- getadmin_gendist_summary(
#   mipanalyzerobject_samples = drcmips$samples,
#   mipanalyzerobject_gendistmat = drcmips$gendistances$IBS_verity_maskhets,
#   type = "mean" # TODO fix matcharg
# )
#
# # Verity IBS force all hets
# drcmips$results$admin_gen_dist$IBS_verity_majallele <- getadmin_gendist_summary(
#   mipanalyzerobject_samples = drcmips$samples,
#   mipanalyzerobject_gendistmat = drcmips$gendistances$IBS_verity_majallele,
#   type = "mean" # TODO fix matcharg
# )
#
# # MalariaGEN DAB
# drcmips$results$admin_gen_dist$DAB_malariagen_wsaf <- getadmin_gendist_summary(
#   mipanalyzerobject_samples = drcmips$samples,
#   mipanalyzerobject_gendistmat = drcmips$gendistances$DAB_malariagen_wsaf,
#   type = "mean" # TODO fix matcharg
# )
#
# # MalariaGEN DAB
# drcmips$results$admin_gen_dist$IBM_0 <- getadmin_gendist_summary(
#   mipanalyzerobject_samples = drcmips$samples,
#   mipanalyzerobject_gendistmat = drcmips$gendistances$IBM_0,
#   type = "mean" # TODO fix matcharg
# )
#
# # IBM05
# drcmips$results$admin_gen_dist$IBM_0 <- getadmin_gendist_summary(
#   mipanalyzerobject_samples = drcmips$samples,
#   mipanalyzerobject_gendistmat = drcmips$gendistances$IBM_0,
#   type = "mean" # TODO fix matcharg
# )
#
# # IBM05
# drcmips$results$admin_gen_dist$IBM_05 <- getadmin_gendist_summary(
#   mipanalyzerobject_samples = drcmips$samples,
#   mipanalyzerobject_gendistmat = drcmips$gendistances$IBM_05,
#   type = "mean" # TODO fix matcharg
# )
#
# # IBD_ghana_drc
# drcmips$results$admin_gen_dist$IBM_05 <- getadmin_gendist_summary(
#   mipanalyzerobject_samples = drcmips$samples,
#   mipanalyzerobject_gendistmat = drcmips$gendistances$MLE_IBD$mle,
#   type = "mean" # TODO fix matcharg
# )
#
## ADMINS
#
# ## make key for liftover
# adminkey_admin1 <- data.frame(
#   item1 = as.character( seq(1, 26) ),
#   admin1_1 = levels(factor(admin_centroids_clst_wght$adm1name))
# )
#
# adminkey_admin2 <- data.frame(
#   item2 = as.character( seq(1, 26) ),
#   admin1_2 = levels(factor(admin_centroids_clst_wght$adm1name))
# )
#
# # Save spatial distance out as a result
# admin_gc_dist.tidy <- broom::tidy(as.dist(drcmips$spatialdist$admin_gc_dist, diag = TRUE,upper=TRUE)) %>%
#   dplyr::mutate(item1 = gsub("admin_long", "", item1),
#                 item2 = gsub("admin_long", "", item2))
#
# admin_gc_dist.tidy <- dplyr::left_join(admin_gc_dist.tidy, adminkey_admin1, by = "item1")
# admin_gc_dist.tidy <- dplyr::left_join(admin_gc_dist.tidy, adminkey_admin2, by = "item2")
#
# # save out
# drcmips$results$admin_gc_dist <- admin_gc_dist.tidy
#
# ## CLUSTER
#
# ## make key for liftover
# clust_key1 <- data.frame(
#   item1 = seq(1,length(clst$hv001)) ,
#   clust_1 = clst$hv001,
#   stringsAsFactors = FALSE
# )
#
# clust_key2 <- data.frame(
#   item2 = seq(1,length(clst$hv001)) ,
#   clust_2 = clst$hv001,
#   stringsAsFactors = FALSE
# )
#
# # Save spatial distance out as a result
# clust_gc_dist.tidy <- broom::tidy(as.dist(drcmips$spatialdist$cluster_gc_dist, diag = TRUE,upper=TRUE))
# clust_gc_dist.tidy <- dplyr::left_join(clust_gc_dist.tidy, clust_key1, by = "item1")
# clust_gc_dist.tidy <- dplyr::left_join(clust_gc_dist.tidy, clust_key2, by = "item2")
#
# # save out
# drcmips$results$cluster_gc_dist <- clust_gc_dist.tidy
#
#
# ###
