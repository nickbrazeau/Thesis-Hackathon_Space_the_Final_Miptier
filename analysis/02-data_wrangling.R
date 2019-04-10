#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle and clean the various covariates in
# the CD2013, weather, and other datasets for the spatial-genetic dist project
#
# Notes: All variables that I will use will either contain a "_fctb/m" or "_cont" to show
#        that I have manipulated/investigated that variable.
#        Men/Women recode combinations (i.e. ha in one and hb in other for same covariate)
#        will be combined to be hab##
#
#        https://dhsprogram.com/pubs/pdf/FR300/FR300.pdf
#
#----------------------------------------------------------------------------------------------------
# libraries and imports
gdrive <- tcltk::tk_choose.dir()

library(tidyverse)
library(osmdata)
source("R/00-base_functions.R")
tol <- 1e-3
set.seed(42)


#--------------------------------------------------------------
# Section 1:Pulling data-map file for all recodes
#--------------------------------------------------------------
bb <- osmdata::getbb("Democratic Republic of the Congo",
                     featuretype = "country",
                     format_out = "sf_polygon")


# cd2013 was under phase 6
# https://dhsprogram.com/publications/publication-DHSG4-DHS-Questionnaires-and-Manuals.cfm
# https://dhsprogram.com/data/Guide-to-DHS-Statistics/ -- note this is version 7
# recode map https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
dt <- readRDS(paste0(gdrive, "/data/raw_data/cd2013_dhs_raw.rds"))

# note subsetting to Mark Janko's "kids" for consistency of IDEEL publications
dt <- dt %>%
  dplyr::filter(markjankosubset == 1)

# drop observations with missing geospatial data
# dt <- dt %>%
#   dplyr::filter(latnum != 0 & longnum != 0) %>%
#   dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>%
#   dplyr::filter(hv103 == 1) # subset to de-facto https://dhsprogram.com/data/Guide-to-DHS-Statistics/Analyzing_DHS_Data.htm

#--------------------------------------------------------------
# Section 2: Looking at recodes, manual data wrangle
#--------------------------------------------------------------

#..........................
# Exposure of Interests
#..........................
# SURVEY CHARACTERISTICS/WEIGHTS
# 1. hv005, cluster level weights
# 2. hv006, Month of interview
# 3. hv007, Year of interview
# 4. hv016, day of interview
#
#
# COINFECTIONS/BIOMARKER VARIABLES
# 1. pfldh coinfection ; (personal + community)
# 2. Hemoglobin Level adjust for altitude and smoking (HA56 & HB56)

# SOCIOECOLOGICAL VARIABLES
# 1. Biological Sex (HV104)
# 2. Age (HV105)
# 3. Main floor material (categorical: HV213)
# 3. Main wall material (categorical: HV214)
# 3. Main roof material (categorical: HV215)
# 3 -> 4. Building material (recode)
# 5. Wealth Index (recode; base wealth is hv270)
# 6. Highest year of education completed (continous: HV108)
# 7. Livestock ownership (categorical: HV246)
# 8. Occupation (categorical: "hv717")
# 9. Number of Household Members (continuous: HV009)
# 10. Number of Children Under Age of 5 (continuous: HV014)


# MALARIA-INTERVENTIONS
# 1. Person slept under an LLIN net (HML20) ; (personal + community)
# 2. Cluster level antimalarial use
#
# PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
# 1. Cluster altitude (HV040)
# 2. Seasonality??
# 3. Temparature
# 4. Precipation
# 6. Type of Residence i.e. Rural or Urban



#########################################################################################################
##############                             SURVEY CHARACTERISTICS                          ##############
#########################################################################################################
#.............
# weights
#.............
dt <- dt %>%
  dplyr::mutate(hv005_wi = hv005/1e6
  )

#.............
# dates
#.............
dt <- dt %>%
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")),
                hvyrmnth_dtmnth = factor(paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")),
                hvyrmnth_dtmnth_lag = factor(paste(lubridate::year(lubridate::rollback(hvdate_dtdmy, roll_to_first = F)),
                                                   lubridate::month(lubridate::rollback(hvdate_dtdmy, roll_to_first = F)), sep = "-"))
  )
xtabs(~dt$hvyrmnth_dtmnth + dt$hvyrmnth_dtmnth_lag)

#.............
# households
#............
summary(dt$hv002) # looks clean if we assume that households are numbered 1 - 34 in each cluster
hs <- dt %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(n = length(hv002),
                   housemax = max(hv002))

summary(hs$housemax) # looks reasonable by cluster

dt <- dt %>%
  dplyr::mutate(houseid = factor(paste0(hv001, "_", hv002)))




#########################################################################################################
##############                          INDIVIDUAL LEVEL VARIABLES                         ##############
#########################################################################################################
#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
#.............
# pfldh/po18s
#.............
summary(dt$pfldh)

dt <- dt %>%
  dplyr::mutate(
    pfldhct_cont = pfctmean,
    pfldhct_cont_log = log(pfldhct_cont + tol),
    pfldhct_cont_log_scale = my.scale(pfldhct_cont_log, center = T, scale = T)
  )

dt[, colnames(dt)[grepl("pfldh", colnames(dt))] ] %>%
  sapply(., summary) # looks clean, all NAs are consistent except Pf but that has to do with double call strategy
# Remember, CT values >cutoff (species dependent) are coded as NA in raw data. Retained here



#.............
# hemoglobin
#.............
levels(factor(haven::as_factor(dt$hc53))) # need to divide by 10

dt <- dt %>%
  dplyr::mutate(hc53_cont = ifelse(hc53 %in% c("992", "993", "996", "999"), NA, hc53),
                hc53_cont = as.numeric(hc53_cont)/10,
                hc53_cont_scale = my.scale(hc53_cont, center = T, scale = T))

summary(dt$hc53_cont)

#..........................................................................................
#                               DEMOGRAPHIC/BEHAVIORAL VARIABLES
#...........................................................................................
#.............
# sex
#.............
levels(factor(haven::as_factor(dt$hv104))) # no missing m/f but still missing factor from haven
sum(is.na(dt$hv104))
dt <- dt %>%
  dplyr::mutate(hv104_fctb = haven::as_factor(hv104),
                hv104_fctb = forcats::fct_recode(hv104_fctb, NULL = "missing"),
                hv104_fctb =  forcats::fct_drop(hv104_fctb),
                hv104_fctb = forcats::fct_rev(forcats::fct_reorder(.f = hv104_fctb, .x = hv104_fctb, .fun = length))
  ) # female to default (b/c 0 and larger SE)

#.............
# age
#.............
levels(factor(haven::as_factor(dt$hc1))) # no missing labels
summary(dt$hc1)
dt <- dt %>%
  dplyr::mutate(hc1_cont = hc1,
                hc1_cont_scale = my.scale(hc1_cont, center = T, scale = T))



#.............
# wealth index
#.............
xtabs(~haven::as_factor(dt$hv270)) # looks OK. Some poor/poorest got placed higher than expected
dt <- dt %>%
  dplyr::mutate(hv270_fctm = haven::as_factor(hv270),
                hv270_fctm = forcats::fct_recode(hv270_fctm, NULL = "missing"),
                hv270_fctm =  forcats::fct_drop(hv270_fctm),
                hv270_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv270_fctm, .x = hv270_fctm, .fun = length))
  ) # poorest to default (b/c sensible and larger SE)


#.............
# Mother's years of education (continuous)
#.............
# NOTE need to recode accounting for hc61 with dhs strat...
# TODO
hist(dt$hc62)
summary(dt$hc62)
dt <- dt %>%
  dplyr::mutate(hc62_cont = ifelse(hc62 %in% c("97", "98", "99"), NA, hc62),
                hc62_cont_scale = my.scale(hc62_cont, center = T, scale = T))

# TODO !!!!!!!



#------------------------------------------
# Family owns livestock, herds, or farm animals
#------------------------------------------
summary(dt$hv246)
table(dt$hv246) # 9 is missing

dt <- dt %>%
  dplyr::mutate(
    hv246_fctb = haven::as_factor(dt$hv246),
    hv246_fctb = forcats::fct_recode(hv246_fctb, NULL = "missing"),
    hv246_fctb =  forcats::fct_drop(hv246_fctb),
    hv246_fctb = forcats::fct_relevel(hv246_fctb, "no")
  )
xtabs(~ dt$hv246 + dt$hv246_fctb, addNA = T)


#------------------------------------------
# children under 5 number
#------------------------------------------
summary(dt$hv014) # looks clean
dt <- dt %>%
  dplyr::mutate(hv014_cont = hv014,
                hv014_cont_scale = my.scale(hv014_cont, center = T, scale = T))

#------------------------------------------
# total household members
#------------------------------------------
summary(dt$hv009) # looks clean
dt <- dt %>%
  dplyr::mutate(hv009_cont = hv009,
                hv009_cont_scale = my.scale(hv009_cont, center = T, scale = T))


#..........................................................................................
#                                 MALARIA-INTERVENTIONS
#..........................................................................................

#.............
# LLIN for INDIVIDUAL
#.............
xtabs(~haven::as_factor(dt$hml10) + haven::as_factor(dt$hml20), addNA = T)
# there are 49 people that slept under ITN" but not LLIN and a few NAs

xtabs(~haven::as_factor(dt$hml19) + haven::as_factor(dt$hml20), addNA = T)
# there are 122 people that slept under "ever treated net" but not LLIN

xtabs(~haven::as_factor(dt$hml10) + haven::as_factor(dt$hml19), addNA = T)
# there are 69 people that slept under "ever treated net" but not a ITN

# BASED on this pattern, am just going to consider LLIN
summary(dt$hml20) # no missing here

dt <- dt %>%
  dplyr::mutate(hml20_fctb = haven::as_factor(hml20))

# #.............
# # LLIN-type of Inseciticide for INDIVIDUAL
# # Note, must have LLIN to have insecticide (120 missing LLIN insecticide types, 8500 no LLIN)
# #.............
# # read insecticide liftover table
# insctcd <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_insecticide_liftover.csv")
#
# dt <- dt %>%
#   dplyr::mutate(hml7 = haven::as_factor(hml7)) %>%
#   left_join(x=., y=insctcd, by="hml7") %>%
#   dplyr::mutate(insctcd_fctm = factor(ifelse(hml20_fctb == "no", "none", insctcd)),
#                 insctcd_fctm = forcats::fct_relevel(insctcd_fctm, "none")
#   )
#
#
# # sanity checks
# xtabs(~dt$insctcd_fctm + dt$hml20_fctb, addNA=T)
# xtabs(~dt$insctcd_fctm + dt$hml7, addNA=T)
# xtabs(~dt$hml20_fctb + dt$hml7, addNA=T)
#
#



#########################################################################################################
##############                             CLUSTER LEVEL VARIABLES                         ##############
#########################################################################################################
# ALL CLUSTER LEVEL VARIABLES WILL BE APPENDED WTIH "_clst"
dtsrvy <- makecd2013survey(survey = dt)

#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
Hbmiss <- dt[is.na(dt$hc53),]
sum(is.na(dt$hc53))
table(Hbmiss$hv001) # looks well dispersed among clusters

coinfxnbiomrk <- dtsrvy %>%
  dplyr::group_by(hv001) %>% # cluster level
  dplyr::summarise(
    pfldh_cont_clst = srvyr::survey_mean(x = pfldh, na.rm = T),
    hc53_cont_clst = srvyr::survey_quantile(hc53_cont, quantiles = c(0.5), vartype = c("se"), na.rm = T)) %>%
  dplyr::mutate(
    pfldh_cont_scale_clst = my.scale(logit(pfldh_cont_clst, tol = tol), center = T, scale = T),
    hc53_cont_scale_clst = my.scale(hc53_cont_clst_q50, center = T, scale = T)
  ) %>%
  dplyr::select(-c(dplyr::ends_with("_se")))

sapply(coinfxnbiomrk, summary) # looks clean, had to remove NAs in Hb because 49 missing (spread across 39 clusters)

dt <- dplyr::left_join(x = dt, y = coinfxnbiomrk)


#..........................................................................................
#                                ECOLOGICAL VARIABLES
#..........................................................................................

#.............
# Precipitation (CHRIPS) & Temperature (Emch/Manny)
#.............
source("R/00-basic_map_functions.R")
# precip data
precip <- list.files(path = paste0(gdrive, "/data/raw_data/weather_data/CHIRPS/", full.names = T))
precipfrst <- lapply(precip, readRasterBB, bb = bb)
precipdf <- tibble::tibble(names = basename(precip)) %>%
  dplyr::mutate(names = gsub("chirps-v2.0.|.tif", "", names),
                year = stringr::str_split_fixed(names, "\\.", n=2)[,1] ,
                mnth =  stringr::str_split_fixed(names, "\\.", n=2)[,2] ,
                hvdate_dtdmy = lubridate::dmy(paste(1, mnth, year, sep = "/")),
                year = lubridate::year(hvdate_dtdmy),
                mnth = lubridate::month(hvdate_dtdmy),
                hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
                precipraster = precipfrst) %>%
  dplyr::select(c("hvyrmnth_dtmnth_lag", "precipraster"))

# temp data
temp <- raster::stack(paste0(gdrive, "data/raw_data/weather_data/emch_manny/cru_ts4.01.2011.2016.tmp.dat.nc")) %>%
  raster::crop(x = ., y = sf::as_Spatial(bb))

tempdf <- tibble::tibble(orignnames = names(temp),
                         names = gsub("X", "", names(temp)),
                         hvdate_dtdmy = lubridate::ymd(names),
                         year = lubridate::year(hvdate_dtdmy),
                         mnth = lubridate::month(hvdate_dtdmy),
                         hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
                         tempraster = lapply(orignnames, raster::subset, x = temp)) %>%
  dplyr::select(c("hvyrmnth_dtmnth_lag", "tempraster"))



wthrnd <- dt[,c("hv001", "hvyrmnth_dtmnth_lag", "geometry", "urban_rura")] %>%
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
wthrnd <- wthrnd[!duplicated(wthrnd),]

wthrnd <- wthrnd %>%
  dplyr::left_join(., tempdf) %>%
  dplyr::left_join(., precipdf)

# Drop in a for loop
wthrnd$precip_lag_cont_clst <- NA
wthrnd$temp_lag_cont_clst <- NA

for(i in 1:nrow(wthrnd)){
  # precip
  wthrnd$precip_lag_cont_clst[i] <-
    raster::extract(x = wthrnd$precipraster[[i]],
                    y = sf::as_Spatial(wthrnd$geometry[i]),
                    buffer = wthrnd$buffer[i],
                    fun = mean,
                    sp = F
    )

  # temp
  wthrnd$temp_lag_cont_clst[i] <-
    raster::extract(x = wthrnd$tempraster[[i]],
                    y = sf::as_Spatial(wthrnd$geometry[i]),
                    buffer = wthrnd$buffer[i],
                    fun = mean,
                    sp = F
    )

}

wthrnd <- wthrnd %>%
  dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag", "precip_lag_cont_clst", "temp_lag_cont_clst")) %>%
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))
sf::st_geometry(wthrnd) <- NULL
dt <- dt %>%
  dplyr::left_join(., wthrnd, by = c("hv001", "hvyrmnth_dtmnth_lag")) %>%
  dplyr::mutate(precip_lag_cont_log_clst = log(precip_lag_cont_clst + tol),
                temp_lag_cont_log_clst = log(temp_lag_cont_clst + tol),
                precip_lag_cont_scale_clst = my.scale(precip_lag_cont_log_clst, center = T, scale = T),
                temp_lag_cont_scale_clst = my.scale(temp_lag_cont_log_clst, center = T, scale = T))



#.............
# Cluster-Level Altitude
#.............
dt <- dt %>%
  dplyr::mutate(alt_dem_cont_clst = ifelse(alt_dem == 9999, NA, alt_dem), # note no missing (likely dropped with missing gps)
                alt_dem_fctb_clst = factor(
                  ifelse(alt_dem_cont_clst > median(alt_dem_cont_clst), "high", "low"),
                  levels = c("low", "high")),
                alt_dem_clst_log = log(alt_dem_cont_clst + tol),
                alt_dem_cont_scale_clst = my.scale(alt_dem_clst_log, center = T, scale = T)
  )


#.............
# Distance to Water Source
#.............
wtrdist_out <- readRDS(paste0(gdrive, "/data/derived_data/hotosm_waterways_dist.rds"))
dt <- dt %>%
  dplyr::left_join(x=., y = wtrdist_out, by = "hv001") %>%
  dplyr::mutate(wtrdist_cont_scale_clst = my.scale(log(wtrdist_cont_clst + tol), center = T, scale = T)
  )


#..........................................................................................
#                           DEMOGRAPHIC/BEHAVIORAL VARIABLES
#..........................................................................................
#.............
# Urbanicity
#.............
xtabs(~haven::as_factor(dt$hv025), addNA = T) # looks OK. Some rural places are now urban, which is consistent with what I expected

dt <- dt %>%
  dplyr::mutate(
    hv025_fctb = haven::as_factor(hv025),
    hv025_fctb = forcats::fct_recode(hv025_fctb, NULL = "missing"), # is none but for safety
    hv025_fctb = forcats::fct_relevel(hv025_fctb, "rural")
  )


#.............
# Distance to Health Site
#.............
hlthdist_out <- readRDS(paste0(gdrive, "/data/derived_data/hotosm_healthsites_dist.rds"))
dt <- dt %>%
  dplyr::left_join(x=., y = hlthdist_out, by = "hv001") %>%
  dplyr::mutate(hlthdist_cont_scale_clst = my.scale(log(hlthdist_cont_clst + tol), center = T, scale = T)
  )

#.............
# LLIN Cluster Usage; median wealth; median educ; health insur
#.............
summary(dt$hml20) # no missing
summary(dt$wlthrcde_fctm)
summary(dt$hc62_cont)

democlust <- dtsrvy %>%
  dplyr::mutate(wlthrcde_fctm_ord = factor(hv270_fctm,
                                           levels = c("poorest", "poorer", "middle", "richer", "richest"),
                                           ordered = T),
                wlthrcde_fctm_ord_num = as.numeric(wlthrcde_fctm_ord)) %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    hml20_cont_clst = srvyr::survey_mean(hml20, vartype = c("se")),
   wlthrcde_fctm_clst = srvyr::survey_median(wlthrcde_fctm_ord_num, quantiles = c(0.5), vartype = c("se")),
    hc62_cont_clst = srvyr::survey_mean(hc62_cont, vartype = c("se"), na.rm = T)
  ) %>%
  dplyr::mutate(
    hml20_cont_scale_clst = my.scale(logit(hml20_cont_clst, tol = tol), center = T, scale = T),
    hc62_cont_scale_clst = my.scale(log(hc62_cont_clst + tol), center = T, scale = T),
    wlthrcde_fctm_scale_clst = factor(floor(wlthrcde_fctm_clst_q50), # taking lower bracket of wealth if median was split
                                         levels = c(1,2,3,4,5),
                                         labels = c("poorest", "poor", "middle", "rich", "richest"))
  ) %>%
  dplyr::select(-c(dplyr::ends_with("_se")))


sapply(democlust, summary) # looks clean

dt <- dplyr::left_join(x = dt, y = democlust)


#.............
# Antimalarial Cluster Usage
#.............
#https://dhsprogram.com/data/Guide-to-DHS-Statistics/
kr <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDKR61FL.rds")
# liftover drug function
# per document, missing goes to NO
missingliftover <- function(x){
  x <- ifelse(x == 9 | is.na(x), 0, x) # 9 is missing
  return(x)
}

denom <- kr %>%
  dplyr::filter(v012 < 60) %>% # less than 5 years
  dplyr::filter(haven::as_factor(b5) == "yes") %>% # currently alive
  dplyr::filter(haven::as_factor(h22) == "yes") %>% # had fever in last two weeks
  dplyr::select(c(paste0("ml13", letters[1:8]), "v001", "v005", "v023")) %>%
  dplyr::select(-c("ml13g")) %>%
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
    fansidar_cont_clst = srvyr::survey_mean(ml13a),
    chloroquine_cont_clst = srvyr::survey_mean(ml13b),
    amodiaquine_cont_clst = srvyr::survey_mean(ml13c),
    quinine_cont_clst = srvyr::survey_mean(ml13d),
    act_cont_clst = srvyr::survey_mean(ml13e),
    otherartm_cont_clst = srvyr::survey_mean(ml13f),
    other_cont_clst = srvyr::survey_mean(ml13h),
    anyatm_cont_clst = srvyr::survey_mean(anyatm)) %>%
  dplyr::mutate(
    fansidar_cont_scale_clst = my.scale(logit(fansidar_cont_clst, tol = tol), center = T, scale = T),
    chloroquine_cont_scale_clst = my.scale(logit(chloroquine_cont_clst, tol = tol), center = T, scale = T),
    amodiaquine_cont_scale_clst = my.scale(logit(amodiaquine_cont_clst, tol = tol), center = T, scale = T),
    quinine_cont_scale_clst = my.scale(logit(quinine_cont_clst, tol = tol), center = T, scale = T),
    act_cont_scale_clst = my.scale(logit(act_cont_clst, tol = tol), center = T, scale = T),
    otherartm_cont_scale_clst = my.scale(logit(otherartm_cont_clst, tol = tol), center = T, scale = T),
    other_cont_scale_clst = my.scale(logit(other_cont_clst, tol = tol), center = T, scale = T),
    anyatm_cont_scale_clst = my.scale(logit(anyatm_cont_clst, tol = tol), center = T, scale = T)
  ) %>%
  dplyr::rename(hv001 = v001)  %>%
  dplyr::mutate(hv001 = as.numeric(hv001)) %>%
  dplyr::select(-c(dplyr::ends_with("_se")))

dt <- dplyr::left_join(dt, kdsrv_fvr_clst, by = "hv001")

dt %>%
  group_by(hv001) %>%
  dplyr::summarise(n = sum(anyatm_cont_clst)) # some clusters missing kids with fevers, so have NAs

#.............
# Cluster Level Kid RDT/Micro
#.............
# https://dhsprogram.com/data/Guide-to-DHS-Statistics/ Percentage of Children Classified in Two Test
pr <- readRDS(paste0(gdrive, "/data/raw_data/dhsdata/datasets/CDPR61FL.rds"))
rdtmicro <- pr %>%
  dplyr::filter(hc1 < 60 & hc1 >= 6) %>% # age
  dplyr::filter(hv042 == 1) %>% # household selected for Hb
  dplyr::filter(hv103 == 1) %>% #de facto
  dplyr::mutate(hv005_wi = hv005/1e6,
                hml32_fctb = haven::as_factor(hml32),
                hml32_fctb = if_else(hml32_fctb %in% c("positive", "negative"),
                                     hml32_fctb, factor(NA)),
                hml32_fctb =  forcats::fct_drop(hml32_fctb),
                hml32_numb = ifelse(hml32_fctb == "positive", 1, ifelse(hml32_fctb == "negative", 0, NA)),

                hml35_fctb = haven::as_factor(hml35),
                hml35_fctb = if_else(hml35_fctb %in% c("positive", "negative"),
                                     hml35_fctb, factor(NA)),
                hml35_fctb =  forcats::fct_drop(hml35_fctb),
                hml35_numb = ifelse(hml35_fctb == "positive", 1, ifelse(hml35_fctb == "negative", 0, NA))
  )


rdtmicro_srvy <- rdtmicro %>% srvyr::as_survey_design(ids = hv001,
                                                      strata = hv023,
                                                      weights = hv005_wi)

rdtmicro_sum <- rdtmicro_srvy %>%
  dplyr::mutate(count = 1) %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    RDTprev_cont_clst = srvyr::survey_mean(hml35_numb, na.rm = T, vartype = c("se")),
    microprev_cont_clst = srvyr::survey_mean(hml32_numb, na.rm = T, vartype = c("se"))
  ) %>%
  dplyr::select(-c(dplyr::ends_with("_se"))) %>%
  dplyr::mutate(RDTprev_cont_scale_clst = my.scale(logit(RDTprev_cont_clst, tol = tol)),
                microprev_cont_scale_clst = my.scale(logit(microprev_cont_clst, tol=tol)),
                hv001 = as.numeric(hv001))

dt <- dplyr::left_join(dt, rdtmicro_sum, by = "hv001")




#..........................................................................................
#                               Final Write Out
#..........................................................................................
saveRDS(dt, file = paste0(gdrive, "/data/derived_data/cd2013_kids_dhs_recode.rds"))


