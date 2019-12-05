#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data and
#  merge PCR results for the Children Datasets.
#  The Pf and the original All Kids data set was made by Mark Janko
#  and are on the Meshnick M Drive
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
# DHS munging
devtools::install_github("OJWatson/rdhs", ref="master")
library(rdhs)


#---------------------------------------------------------------------------------
# Using rDHS to pull down CD2013
# This will be used for the both kids and the adults
#---------------------------------------------------------------------------------
# https://ojwatson.github.io/rdhs/
rdhs::set_rdhs_config(email = "nbrazeau@med.unc.edu",
                      project = "Malaria Spatiotemporal Analysis",
                      config_path = "rdhs.json",
                      global = FALSE,
                      cache_path = "~/Documents/GitHub/Space_the_Final_Miptier/data/raw_data/dhsdata/")

survs <- rdhs::dhs_surveys(countryIds = c("CD"),
                           surveyYearStart = 2013)


datasets <- rdhs::dhs_datasets(surveyIds = survs$SurveyId,
                               fileFormat = "flat")

# download all DHS datsets
downloads <- rdhs::get_datasets(datasets$FileName)


#.........................................
# Read and Set up DHS Data
#.........................................
pr <- readRDS(file = "data/raw_data/dhsdata/datasets/CDPR61FL.rds")

# fix kids barcode
pr <- pr %>%
  mutate(hmid = paste(hv001, hv002, hvidx, sep = "_"),
         barcode = gsub(" ", "", sh312),
         barcode = tolower(barcode),
         barcode = ifelse(barcode %in% c("99993", "99994", "99995", "99996"), NA, barcode),
         barcode = ifelse(grepl("\\?", barcode), NA, barcode)
  )



#---------------------------------------------------------------------------------
# Pull in Kids qPCR Results
#---------------------------------------------------------------------------------
pfpvpcr.kids <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Children_Construction/PfPv_allchildren_v1.csv",
                                col_names = T) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(barcode = sh312,
                pvctcrrct = ct_corrected)


# note, these 7,250 kids are the under-5s versus the over-5s, match Mark's numbers perfectly
xtabs(~pfpvpcr.kids$pr_present, addNA = T)

# going to subset to note "missing data" children versus
# "full data" children
# and bring in PR
dt <- pfpvpcr.kids %>%
  dplyr::mutate(U5_O5 = ifelse(!is.na(pr_present), "full", "missing")) %>%
  dplyr::left_join(x=., y=pr, by = "barcode")




#---------------------------------------------------------------------------------
# DHS spatial shapes/data
#---------------------------------------------------------------------------------

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_KA_Compare/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>%
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR


dt <- left_join(dt, ge, by = "hv001")

# drop observations with missing geospatial data
dt.geo <- dt %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

spacemipsge <- dt.geo %>%
  dplyr::select(c("adm1name", "hv001", "latnum", "longnum", "geometry")) %>%
  dplyr::filter(!duplicated(.))

# save out the final product
saveRDS(spacemipsge, file = "data/derived_data/spacemips_GE.rds")
saveRDS(dt.geo, file = "data/derived_data/DHS_qPCR_allkids_geo.rds")



















