# per cluster data frame

space <- readRDS(paste0(gdrive,"/data/derived_data/cd2013_gen_space_epi_final.rds"))
meta <- readRDS(paste0(gdrive,"/data/derived_data/cd2013_kids_dhs_recode.rds"))
gen <- readRDS(paste0(gdrive,"/data/derived_data/cd2013_dhs_bigbarcodemips.rds"))

gen$samples$ADM1NAME <- factor(as.character(gen$samples$ADM1NAME))
gen$samples$DHSREGNA <- factor(as.character(gen$samples$DHSREGNA))

# grab the eco variables as needed
gen$samples <- left_join(gen$samples,
                         meta[,c("hv001","hv005","hc1","dhsregna","temp_lag_cont_clst",
                                 "annual_precipitation_2015","alt_dem_cont_clst","travel_times_2015",
                                 "anyatm_cont_clst","microprev_cont_clst","latnum","longnum",
                                 "RDTprev_cont_clst","pfldh_cont_clst")],
                         by = "hv001")

# group down
clust_meta <- group_by(gen$samples, hv001) %>%
  summarise(admin1 = unique(dhsregna),
            rdt_prev = mean(RDTprev_cont_clst),
            micro_prev = mean(microprev_cont_clst),
            ldh_prev = mean(pfldh_cont_clst),
            drug_use = mean(anyatm_cont_clst),
            precip = mean(annual_precipitation_2015),
            altitude = mean(alt_dem_cont_clst),
            travel = mean(travel_times_2015),
            lat = mean(latnum),
            long = mean(longnum)
            )


# previously calculated admin survey weighted variables
meta_admin <- readRDS(paste0(gdrive,"/data/derived_data/cd2013_kids_dhs_admin1_recode.rds"))

rm(gen); gc()
gen <- readRDS(paste0(gdrive,"/data/derived_data/cd2013_dhs_bigbarcodemips.rds"))

gen$samples$ADM1NAME <- factor(as.character(gen$samples$ADM1NAME))
gen$samples$DHSREGNA <- factor(as.character(gen$samples$DHSREGNA))


gen$samples <- left_join(gen$samples,
                         meta_admin[,c("hv001","hv005","hc1","dhsregna","temp_lag_cont_clst",
                                 "annual_precipitation_2015","alt_dem_cont_clst","travel_times_2015",
                                 "anyatm_cont_clst","microprev_cont_clst",
                                 "RDTprev_cont_clst","pfldhct_cont")],
                         by = "hv001")

admin_meta <- group_by(gen$samples, ADM1NAME) %>%
  summarise(
    dhs_region = unique(dhsregna),
    micro_prev = mean(microprev_cont_clst, na.rm=TRUE),
    rdt_prev = mean(RDTprev_cont_clst, na.rm=TRUE),
    ldh_prev = mean(pfldhct_cont, na.rm=TRUE),
    drug_use = mean(anyatm_cont_clst, na.rm=TRUE),
    precip = mean(annual_precipitation_2015, na.rm=TRUE),
    altitude = mean(alt_dem_cont_clst, na.rm=TRUE),
    travel = mean(travel_times_2015, na.rm=TRUE),
    lat = mean(LATNUM),
    long = mean(LONGNUM)
  )
saveRDS(admin_meta, paste0(gdrive,"/data/derived_data/cd2013_admin_covariates.rds"))
