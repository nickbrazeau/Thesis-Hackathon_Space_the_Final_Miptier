#................................................................................................................................
# Purpose of this script is to pull together covariates
#................................................................................................................................
library(tidyverse)
library(sf)
source("R/basics.R")

# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>%
  dplyr::rename(hv001 = dhsclust)

#.............
# Pf Prevalence
#.............
rdtliftover <- readr::read_csv("internal_datamap_files/RDT_liftover.csv")
microliftover <- readr::read_csv("internal_datamap_files/microscopy_liftover.csv")
pr.prev <- readRDS("data/derived_data/DHS_qPCR_allkids_geo.rds") %>%
  dplyr::select(c("hv001", "hml32", "hml35", "pfldh"))
pr.prev <- pr.prev %>%
  dplyr::left_join(., y=rdtliftover, by = "hml35") %>%
  dplyr::left_join(., y=microliftover, by = "hml32") %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    n = n(),
    qPCRprev = mean(pfldh, na.rm = T),
    qPCRse = sd(pfldh, na.rm = T)/sqrt(n),

    rdtprev = mean(RDTresult, na.rm = T),
    rdtse = sd(RDTresult, na.rm = T)/sqrt(n),

    microprev = mean(microscopyresult, na.rm = T),
    microprevse = sd(microscopyresult, na.rm = T)/sqrt(n)

  )

pr.prev.scaled <- pr.prev %>%
  dplyr::mutate(
    qPCRprev_scale = my.scale(logit(qPCRprev, tol = 1e-3)),
    rdtprev_scale =  my.scale(logit(rdtprev, tol = 1e-3)),
    microprev_scale =  my.scale(logit(microprev, tol = 1e-3))
  ) %>%
  dplyr::select(c("hv001", "qPCRprev_scale", "rdtprev_scale", "microprev_scale")) %>%
  dplyr::mutate(hv001 = as.numeric(hv001))



#.............
# Net Use
#.............
netuse <- readRDS("data/derived_data/DHS_kids_net_use.rds") %>%
  dplyr::select(c("hv001", "ITNperc_logit_scale"))

#.............
# ACT Use
#.............
rxuse <- readRDS("data/derived_data/DHS_kids_act_use_imputed.rds") %>%
  dplyr::select("hv001", "anyatmperc_logit_scale")

#.............
# Distance Hospital
#.............
hospdist <- readRDS("data/derived_data/hlthdist_out_minduration.rds") %>%
  dplyr::mutate(hlthdist_log_scale = my.scale(log(hlthst_nrst_duration))) %>%
  dplyr::select(-c("hlthst_nrst_duration"))

#.............
# Urbanicity
#.............
urb <- readRDS("data/derived_data/DRC_urbanicity.rds")

#.............
# Weather
#.............
wthr <- readRDS("data/derived_data/vividep_weather_recoded_mean.rds") %>%
  dplyr::mutate(precip_mean_scale = my.scale(precip_mean_cont_clst),
                temp_mean_scale = my.scale(temp_mean_cont_clst)) %>%
  dplyr::select(-c("precip_mean_cont_clst", "temp_mean_cont_clst"))

#.............
# Cluster-Level Altitude
#.............
alt <- ge %>%
  dplyr::select("hv001", "alt_dem") %>%
  dplyr::mutate(alt_dem_log_scale = my.scale(log(alt_dem))) %>%
  dplyr::select(-c("alt_dem"))


#..........................................................
# Join and save
#..........................................................
covariatematrix <- dplyr::left_join(pr.prev.scaled, netuse, by ="hv001") %>%
  dplyr::left_join(., rxuse, by = "hv001") %>%
  dplyr::left_join(., hospdist, by = "hv001") %>%
  dplyr::left_join(., urb, by = "hv001") %>%
  dplyr::left_join(., wthr, by = "hv001") %>%
  dplyr::left_join(., alt, by = "hv001")

# 492 to 490 because clusters 70 and 510 didn't have under 5s?


saveRDS(object = covariatematrix, "data/derived_data/covariatematrix.rds")


