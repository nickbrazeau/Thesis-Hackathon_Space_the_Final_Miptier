#................................................................
# Purpose of this script is to calculate malaria intervention
# by cluster in the DRC
#................................................................
library(tidyverse)
library(sf)
source("R/basics.R")
tol <- 1e-3

# Notes on Lonely PSUs
# http://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
options(survey.lonely.psu="adjust")



#...........................................................
# Bednet Cluster Usage
#...........................................................
#https://dhsprogram.com/data/Guide-to-DHS-Statistics/
# Use of Mosquito Nets by Children
pr <- readRDS(file = "data/raw_data/dhsdata/datasets/CDPR61FL.rds")
pr.defacto.U5 <- pr %>%
  dplyr::mutate(hv005 = hv005/1e6,
                ITN = ifelse(hml16 %in% c(1,2), 1, 0), # "Nets with missing information on type are considered not to be ITNs."
                ITN_fctb = factor(ITN, levels = c(0,1), labels = c("No", "Yes"))) %>%
  dplyr::filter(hv103 == 1 & hml16 < 5)


kdsrv_srvyr <- pr.defacto.U5 %>% srvyr:::as_survey_design(ids = hv001,
                                                          strata = hv023,
                                                          weights = hv005)

kdsrv_nets <- kdsrv_srvyr %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    ITNperc = srvyr::survey_mean(ITN)
    )

xtabs(~pr.defacto.U5$ITN_fctb + pr.defacto.U5$ITN) # looks clean


# NOW SCALE
kdsrv_nets <- kdsrv_nets %>%
  dplyr::mutate(
    ITNperc_logit = logit(ITNperc, tol = tol),
    ITNperc_logit_scale = my.scale(ITNperc_logit, center = T, scale = T)
  )


saveRDS(file = "data/derived_data/DHS_kids_net_use.rds", object = kdsrv_nets)



#...........................................................
# Antimalarial Cluster Usage
#...........................................................
#https://dhsprogram.com/data/Guide-to-DHS-Statistics/
kr <- readRDS(file = "data/raw_data/dhsdata/datasets/CDKR61FL.rds")

# liftover drug function
# per document, missing goes to NO
missingliftover <- function(x){
  x <- ifelse(x == 9 | is.na(x), 0, x) # 9 is missing
  return(x)
}


denom <- kr %>%
  dplyr::filter(b8 < 5) %>% # less than 5 years -- dhs calls for b19<60, but no b19 variable in kr, ir, hr... this should be sufficient as it is less than 5 years
  dplyr::filter(haven::as_factor(b5) == "yes") %>% # currently alive
  dplyr::filter(haven::as_factor(h22) == "yes") %>% # had fever in last two weeks
  dplyr::select(c(paste0("ml13", letters[1:8]), "v001", "v005", "v023")) %>%
  dplyr::select(-c("ml13g")) %>% # coded as NA in cd2013
  dplyr::mutate(v005 = v005/1e6)

# clean up
denom[, grepl("ml13", colnames(denom))] <- lapply(denom[, grepl("ml13", colnames(denom))],
                                                  missingliftover) %>% dplyr::bind_cols(.)
# check
sapply(denom, summary)
# add any use in
denom$anyatm = as.numeric( apply(denom[,grepl("ml13", colnames(denom))],
                                 1, function(x){return(any( x == 1))}) )
# Note, some individuals took multiple drugs. OK because small percent 36/1560

kdsrv_fvr <- denom %>% srvyr:::as_survey_design(ids = v001,
                                                strata = v023,
                                                weights = v005)

kdsrv_fvr_clst <- kdsrv_fvr  %>%
  dplyr::mutate(count = 1) %>%
  dplyr::group_by(v001) %>%
  dplyr::summarise(
    n = srvyr::survey_total(count),
  anyatmperc = srvyr::survey_mean(anyatm)) %>%
  dplyr::rename(hv001 = v001)  %>%
  dplyr::mutate(hv001 = as.numeric(hv001)) %>%
  dplyr::select(-c(dplyr::ends_with("_se")))



#-----------------------------------------------------------
# Not all clusters had kids with fever in the last 2 weeks,
# Need to impute 8 missing clusters
#-----------------------------------------------------------
# find missing clusters
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::select(c("hv001", "latnum", "longnum", "geometry"))


kdsrv_fvr_clst.imp <- dplyr::left_join(ge, kdsrv_fvr_clst, by = "hv001")

sum(is.na(kdsrv_fvr_clst.imp$anyatmperc))

#.................
# Look to see if KNN is reasonable
#.................
mssngclust <- kdsrv_fvr_clst.imp %>%
  dplyr::filter(is.na(n)) %>%
  dplyr::select(-c("latnum", "longnum"))

notmssngclust <- kdsrv_fvr_clst.imp %>%
  dplyr::filter(!is.na(n)) %>%
  dplyr::select(-c("latnum", "longnum"))

ggplot() +
  geom_sf(data=notmssngclust, color= "black") +
  geom_sf(data=mssngclust, color= "red", size = 2, alpha  = 0.5)

# going to have border issues but that's OK

find_mean_actuse_five_clusters <- function(missingcluster.sf, knownclusters.sf){
  # obviously this is not a general function to be exported
  dist <- raster::pointDistance(p1 = sf::as_Spatial(missingcluster.sf),
                                p2 = sf::as_Spatial(knownclusters.sf),
                                lonlat = T)

  # find 5 nearby clusters
  dist.sorted.5 <- sort(dist)[1:5]
  nrbyclstrs <- which(dist %in% dist.sorted.5)
  # sanity check
  if(length(nrbyclstrs) != 5){
    stop("nearest clusters don't match up")
  }
  nrbyclstrs <- knownclusters.sf$hv001[nrbyclstrs]

  ret <- knownclusters.sf %>%
    dplyr::filter(hv001 %in% nrbyclstrs) %>%
    dplyr::summarise(
      anyatmperc = mean(anyatmperc)
    )

  sf::st_geometry(ret) <- NULL # don't need spatial points anymore

  missingcluster.sf$anyatmperc <- ret$anyatmperc

  sf::st_geometry(missingcluster.sf) <- NULL # don't need spatial points anymore...
  return(missingcluster.sf)
}


# run it all
mssngclust.list <- mssngclust %>%
  base::split(1:nrow(.))


mssngclust <- lapply(mssngclust.list, find_mean_actuse_five_clusters, knownclusters.sf = notmssngclust) %>%
  dplyr::bind_rows(.)


# Combine to useful parts
sf::st_geometry(notmssngclust) <- NULL # don't need spatial points anymore :)

kdsrv_fvr_clst.imp <- dplyr::bind_rows(notmssngclust, mssngclust) %>%
  dplyr::arrange(hv001)



# NOW SCALE
kdsrv_fvr_clst.imp <- kdsrv_fvr_clst.imp %>%
  dplyr::mutate(
    anyatmperc_logit = logit(anyatmperc, tol = tol),
    anyatmperc_logit_scale = my.scale(anyatmperc_logit, center = T, scale = T)
  )  %>%
  dplyr::select(-c("n"))



saveRDS(file = "data/derived_data/DHS_kids_act_use_imputed.rds", object = kdsrv_fvr_clst.imp)




