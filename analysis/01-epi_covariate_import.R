#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data and merge DHS recodes for CD2013 and
# various other covariates needed for our space-genetic-distance analysis
# Will then merge this to the hlth-clinic/water distance and climate data in the 01-hotosm_chirps file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
# DHS munging
devtools::install_github("OJWatson/rdhs", ref="master")
library(rdhs)

gdrive <- tcltk::tk_choose.dir()

#---------------------------------------------------------------------------------
# Using rDHS to pull down CD2013
#---------------------------------------------------------------------------------
# https://ojwatson.github.io/rdhs/
rdhs::set_rdhs_config(email = "nbrazeau@med.unc.edu",
                      project = "Malaria Spatiotemporal Analysis",
                      config_path = "rdhs.json",
                      global = FALSE,
                      cache_path = paste0(gdrive, "/data/raw_data/dhsdata/"))

survs <- rdhs::dhs_surveys(countryIds = c("CD"),
                           surveyYearStart = 2013)


datasets <- rdhs::dhs_datasets(surveyIds = survs$SurveyId,
                               fileFormat = "flat")




# download all DHS datsets
downloads <- rdhs::get_datasets(datasets$FileName)


#---------------------------------------------------------------------------------
# Read in and match PR barcode to qpcr data
#---------------------------------------------------------------------------------
pr <- readRDS(file = paste0(gdrive, "/data/raw_data/dhsdata/datasets/CDPR61FL.rds"))
# note can circle back and use the MR and IR recodes for other cluster level data later if we want

#---------------------------------------------------------------------------------
# DHS spatial shapes/data
#---------------------------------------------------------------------------------

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = paste0(gdrive, "/data/raw_data/dhsdata/datasets/CDGE61FL.rds")))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>%
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR

gepr <- left_join(x=pr, y=ge, by = "hv001")

#---------------------------------------------------------------------------------
# DHS Geospatial Covariates
#---------------------------------------------------------------------------------
gc <- readRDS(paste0(gdrive, "/data/raw_data/dhsdata/datasets/CDGC62FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR

gcgepr <- dplyr::left_join(gepr, gc, by = c("dhsid", "dhscc", "dhsyear", "hv001"))



#---------------------------------------------------------------------------------
# FINAL -- subset the DHS data to just the PCR data
#---------------------------------------------------------------------------------

# Mark created this file for the Janko et al. 2017 MS and is stored on the M-drive
pfpcr <- readr::read_csv(file="/Volumes/share/1. Data/3. Data Sets/2013-14 Kids Database/MappingDataFullKids.csv",
                         col_names = T) %>%
  dplyr::select(c("sh312", "result", "pf_meanCt")) %>%
  dplyr::rename(pfldh = result,
                pfctmean = pf_meanCt) %>%
  dplyr::mutate(sh312 = tolower(sh312),
                sh312 = gsub(" ", "", sh312),
                markjankosubset = 1)


gcgepr <- gcgepr %>%
  dplyr::mutate(sh312 = tolower(sh312),
                sh312 = gsub(" ", "", sh312)) %>%
  dplyr::left_join(x = ., y = pfpcr, by = "sh312")

# sanity
snty <- dplyr::inner_join(gcgepr, pfpcr, by  = "sh312")
if(nrow(snty) != sum(!is.na(gcgepr$pfldh))){
  stop("There are missing dhs barcodes in the kids full data")
}


#---------------------------------------------------------------------------------
# write out
#---------------------------------------------------------------------------------
# write out joined PR, IR, MR, GE, and GC recodes
# in order to perserve sf features, covert to sf (and dataframe)
gcgepr <- sf::st_as_sf(gcgepr)



saveRDS(gcgepr, file = paste0(gdrive, "/data/raw_data/cd2013_dhs_raw.rds"))


