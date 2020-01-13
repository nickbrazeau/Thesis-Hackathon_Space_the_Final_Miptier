#..............................................................
# Purpose of this script is to merge in the PR recode with
# the samples that we actually sequenced
#..............................................................
library(tidyverse)
#..............................................................
# Read In Sequenced Samples
#..............................................................
mtdt <- readRDS("data/derived_data/sample_metadata.rds")

#..............................................................
# Read PR
#..............................................................
dt.geo <- readRDS("data/derived_data/DHS_qPCR_DRCsmpls_geo.rds")

#..............................................................
# Subset to Sequenced samples
#..............................................................
dt <- dt.geo %>%
  dplyr::filter(barcode %in% tolower(mtdt$barcode))

########################################################################
########################################################################
################          WRANLGE IND COVARS            ################
#########################################################################
########################################################################

#..........................
# Exposure of Interests
#..........................
# SURVEY CHARACTERISTICS/WEIGHTS
# 2. hv002, households

# SOCIOECOLOGICAL VARIABLES
# 1. Biological Sex (HV104)
# 2. Age (HV105)
# 3. Main floor material (categorical: HV213)
# 3. Main wall material (categorical: HV214)
# 3. Main roof material (categorical: HV215)
# 3 -> 4. Building material (recode)
# 5. Wealth Index (my recode; base wealth is hv270)
# 6. Highest year of education completed (continous: HV108)
# 7. Occupation (categorical: "hv717")
# 8. Number of Household Members (continuous: HV009)


# MALARIA-INTERVENTIONS
# 1. ITN Use
#

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


#.............
# dates
#.............
dt <- dt %>%
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")),
                hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-"))



#########################################################################################################
##############                          INDIVIDUAL LEVEL VARIABLES                         ##############
#########################################################################################################
#.............
# sex
#.............
levels(factor(haven::as_factor(dt$hv104))) # no missing m/f but still missing factor from haven
sum(is.na(dt$hv104))
dt <- dt %>%
  dplyr::mutate(sex_fctb = haven::as_factor(hv104),
                sex_fctb = forcats::fct_recode(sex_fctb, NULL = "missing"),
                sex_fctb =  forcats::fct_drop(sex_fctb),
                sex_fctb = forcats::fct_rev(forcats::fct_reorder(.f = sex_fctb, .x = sex_fctb, .fun = length))
  ) # female to default (b/c 0 and larger SE)


#.............
# age
#.............
dt <- dt %>%
  dplyr::mutate(age = ifelse(U5_O5 == "adult" & haven::as_factor(hv104) == "female",
                                   ha1, hb1),
                age = ifelse(U5_O5 == "full", hc1/12, age),
                age = ifelse(U5_O5 == "missing", 5.1, age), # impute these kids ages at 5.1
                age = ifelse(age %in% c("97", "98", "99"), NA, age)
                )

summary(dt$age)
hist(dt$age)

#.............
# main floor, wall, roof
#.............
# recode to rudimentary or non-rudimentary per previous manuscript (PMID: 28222094)
# then recode house to modern or traditional (final covar)
# https://pdfs.semanticscholar.org/e290/cf81bdb182696505952f37d1c910db86925a.pdf

# floor
floor <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_floor_liftover.csv")
dt <- dt %>%
  dplyr::mutate(hv213 = haven::as_factor(hv213),
                hv213 =  forcats::fct_drop(hv213))
dt <- dt %>%
  left_join(x=., y=floor, by="hv213") %>%
  dplyr::mutate( hv213_liftover = factor(floortype) )

xtabs(~dt$hv213 + dt$hv213_liftover, addNA = T)


# wall
wall <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_wall_liftover.csv")
dt <- dt %>%
  dplyr::mutate(hv214 = haven::as_factor(hv214),
                hv214 =  forcats::fct_drop(hv214))
dt <- dt %>%
  left_join(x=., y=wall, by="hv214") %>%
  dplyr::mutate( hv214_liftover = factor(walltype) )

xtabs(~dt$hv214 + dt$hv214_liftover, addNA = T)


# roof
roof <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_roof_liftover.csv")
dt <- dt %>%
  dplyr::mutate(hv215 = haven::as_factor(hv215),
                hv215 =  forcats::fct_drop(hv215))
dt <- dt %>%
  left_join(x=., y=roof, by="hv215") %>%
  dplyr::mutate( hv215_liftover = factor(rooftype) )

xtabs(~dt$hv215 + dt$hv215_liftover, addNA = T)

# Final Liftover to a binary variable
# 6 missing in floor; 5 missing in walls
wallfloormiss <- dt %>%
  dplyr::select(c(hv213, hv214)) %>%
  dplyr::mutate(hv213 = haven::as_factor(dt$hv213),
                hv214 = haven::as_factor(dt$hv214)) %>%
  dplyr::filter(hv213 == "missing" | hv214 == "missing")
xtabs(~wallfloormiss$hv213 + wallfloormiss$hv214, addNA = T)
# different observations missing for floor v. wall


# Based on recent evidence, I think all metal roofs should be considered modern
# because of the indoor temperature causing mosquito death PMC6533302
wallfloorroofrecode <- dt %>%
  dplyr::select(c("hv213_liftover", "hv214_liftover", "hv215_liftover"))
sf::st_geometry(wallfloorroofrecode) <- NULL
dt <- dt %>%
  dplyr::mutate(
    housecount = apply(wallfloorroofrecode, 1, function(x){return(sum(x == "non-rudimentary"))}),
    hv21345_fctb = ifelse(housecount == 3, "modern", "traditional"), # per PMID: 28222094
    hv21345_fctb = ifelse(haven::as_factor(dt$hv215) == "metal", "modern", hv21345_fctb), # per PMC6533302
    hv21345_fctb = factor(hv21345_fctb),
    hv21345_fctb = relevel(hv21345_fctb, "modern"))


# check -- seems reasonable
xtabs(~dt$hv21345_fctb, addNA = T)

dt <- dt %>%
  dplyr::rename(housetype = hv21345_fctb)


#.............
# wealth index
#.............
dt <- dt %>%
  dplyr::rename(wealth = hv270)


#..........................................................................................
#                                 MALARIA-INTERVENTIONS
#..........................................................................................

#.............
# ITN for INDIVIDUAL
#.............
# https://dhsprogram.com/Data/Guide-to-DHS-Statistics/index.cfm
# Use of Mosquito Nets by Persons in the Household
# Going to use the definition by Tusting et al. 2017 PMC5319641

xtabs(~haven::as_factor(dt$hml12)) # note, we have no one who slept under a both treated and untreated net
# and based on above, we have no missing net information

xtabs(~haven::as_factor(dt$hml12))
xtabs(~haven::as_factor(dt$hml10))

# using PMC5319641 definition
dt <- dt %>%
  dplyr::mutate(
    ITN_fctb = ifelse(
      # (i) long-lasting insecticidal nets that were <= 3 y old at the time of survey
      !(haven::as_factor(hml4) %in% c("more than 3 years ago", "don't know", "missing")) & haven::as_factor(hml20) == "yes", 1,
      # (ii) conventional ITNs that were 1 y old
      ifelse(haven::as_factor(hml4) %in% c(0:12) & haven::as_factor(hml12) == "only treated (itn) nets", 1,
             # or were retreated within the year before the survey
             ifelse(haven::as_factor(hml9) %in% c(0:12) & haven::as_factor(hml12) == "only treated (itn) nets", 1,
                    0))), # we know no missing net from above. they either reported yes or no at some level
    ITN_fctb = factor(ITN_fctb, levels = c(0,1), labels = c("no", "yes")),
    ITN_fctb = forcats::fct_drop(ITN_fctb),
    ITN_fctb = relevel(ITN_fctb, "yes")
  )


# sanity check
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml10))
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml12))
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml20))


#.............
# LLIN-Net
#.............
dt <- dt %>%
  dplyr::mutate(hml20_fctb = haven::as_factor(hml20)) %>%
  dplyr::rename(LLIN_fctb = hml20_fctb)


#..............................................................
# Save out
#..............................................................
dt <- dt %>%
  dplyr::rename(ind_survey_lvl = U5_O5)
indvars <- c("hv001", "barcode", "ind_survey_lvl")
covars <- c("houseid", "hvdate_dtdmy", "hvyrmnth_dtmnth", "sex_fctb", "age", "housetype",
            "wealth", "ITN_fctb", "LLIN_fctb")
dt.ind <- dt %>%
  dplyr::select(c(indvars, covars))
saveRDS(dt.ind, file = "data/derived_data/dhs_ind_covars.RDS")




