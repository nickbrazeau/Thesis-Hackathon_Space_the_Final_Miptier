#..............................................................
# Purpose of this script is to find DRC samples and assoc
# metadata
#..............................................................

bigbarcode <- readRDS(file = "data/raw_data/bigbarcode_genetic_data/biallelic_processed.rds")
smplmtdt <- bigbarcode$samples
smplmtdt <- smplmtdt %>%
  dplyr::filter(Country == "DRC") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(name = id)
saveRDS(smplmtdt, file = "data/derived_data/sample_metadata.rds")



