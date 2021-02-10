## .................................................................................
## Purpose: Look at COI estimation from RMCL from big barcode manuscript
##
## Notes:
## .................................................................................
library(tidyverse)
#............................................................
##### Data Wrangle ####
#...........................................................
coi <- readRDS("data/raw_data/RMCL_results/non_summariased_cois.rds") %>%
  dplyr::filter(region_denom == 1) %>% # just subset to DRc
  dplyr::filter(region == "DRC") %>%
  dplyr::filter(gt == 0.1) %>% # same as in Big Barcode Manuscript
  dplyr::mutate(barcode = ifelse(nchar(name) == 9, paste(strsplit(name, split = "")[[1]][5:nchar(name)], collapse = ""),
                                 ifelse(nchar(name) == 8, paste(strsplit(name, split = "")[[1]][4:nchar(name)], collapse = ""),
                                        ifelse(nchar(name) == 6, paste(strsplit(name, split = "")[[1]][2:nchar(name)], collapse = ""),
                                               name))))


# subset
monoclonals <- coi %>% dplyr::filter(median == 1)
dir.create("data/derived_data/RMCL_results/", recursive = T)
saveRDS(monoclonals, "data/derived_data/RMCL_results/monoclonals.rds")

