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

dir.create("data/raw_data/dhsdata", recursive = T)
#---------------------------------------------------------------------------------
# Using rDHS to pull down CD2013
# This will be used for the both kids and the adults
#---------------------------------------------------------------------------------
rdhs::set_rdhs_config(email = "nbrazeau@med.unc.edu",
                      project = "Malaria Spatiotemporal Analysis",
                      config_path = "rdhs.json",
                      global = FALSE,
                      cache_path = "data/raw_data/dhsdata/")

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
dt.kids <- pfpvpcr.kids %>%
  dplyr::mutate(U5_O5 = ifelse(!is.na(pr_present), "full", "missing")) %>%
  dplyr::left_join(x=., y=pr, by = "barcode")

#..............................................................
# missing children cluster information
#..............................................................
missclustspace <- readRDS("data/raw_data/missing_children_cluster_info/clean_meta_fromOJ.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("barcode", "hv001")) %>%
  dplyr::rename(missingkidhv001 = hv001) %>%
  dplyr::mutate(barcode = tolower(barcode))

dt.kids <- dt.kids %>%
  dplyr::left_join(., missclustspace, by = "barcode") %>%
  dplyr::mutate(
    hv001 = ifelse(U5_O5 == "missing", missingkidhv001, hv001)
  )


#..............................................................
# Pull in Adult Results and Barcodes
#..............................................................
# read in data
pfpcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pf_alladults_v4.csv",
                         col_names = T) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("hivrecode_barcode", "pfldh", "fcq_mean")) %>%
  dplyr::rename(pfctmean = fcq_mean)

ar <- readRDS(file = "data/raw_data/dhsdata/datasets/CDAR61FL.rds") %>%
  dplyr::rename(hv001 = hivclust,
                hv002 = hivnumb,
                hvidx = hivline,
                hivrecode_barcode = hiv01) %>%
  dplyr::mutate(hivrecode_barcode = gsub(" ", "", hivrecode_barcode),
                hivrecode_barcode = tolower(hivrecode_barcode)) # rename and fix barcode for ar
# match HIV/PCR barcodes with the PR recode
arpr <- dplyr::inner_join(ar,pr, by = c("hv001", "hv002", "hvidx"))
arpr <- arpr %>%
  dplyr::mutate(hmid = paste(hv001, hv002, hvidx, sep = "_")) %>%
  dplyr::select(-c("barcode")) # drop creation of PR barcode above because we need AR barcode now

dt.adults <- arpr %>%
  dplyr::left_join(pfpcr, ., by = "hivrecode_barcode") %>%
  dplyr::rename(barcode = hivrecode_barcode) %>%
  dplyr::mutate(U5_O5 = "adult") # for missing children tracking

#..............................................................
# Merge Kids and Adults
#..............................................................
drop.adults <- colnames(dt.adults)[!colnames(dt.adults) %in% colnames(dt.kids)] # adult stuff relates to HIV
drop.kids <- colnames(dt.kids)[!colnames(dt.kids) %in% colnames(dt.adults)] # kid stuff relates to vivax

dt.adults <- dt.adults[,c( !colnames(dt.adults) %in% c(drop.adults, drop.kids)) ]
dt.kids <- dt.kids[,c( !colnames(dt.kids) %in% c(drop.adults, drop.kids)) ]
# bind
dt <- rbind.data.frame(dt.adults, dt.kids)

#---------------------------------------------------------------------------------
# DHS spatial shapes/data
#---------------------------------------------------------------------------------

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
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

#..............................................................
# Save out
#..............................................................
dir.create("data/derived_data", recursive = T)
saveRDS(spacemipsge, file = "data/derived_data/spacemips_GE.rds")
saveRDS(dt.geo, file = "data/derived_data/DHS_qPCR_DRCsmpls_geo.rds")


