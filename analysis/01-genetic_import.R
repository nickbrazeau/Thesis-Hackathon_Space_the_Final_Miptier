#----------------------------------------------------------------------------------------------------
# Purpose of this script is simply to filter the genetic data that was previously
# filtered and included in the Big Barcode paper by Bob Verity
#
#----------------------------------------------------------------------------------------------------


# libraries and imports
gdrive <- tcltk::tk_choose.dir()
library(tidyverse)


# Bob Verity filtered this for the
# big barcode paper
mipbigpanel <- readRDS(paste0(gdrive, "/data/raw_data/mip_genetic_data_from_Bob_big_barcode/biallelic_processed.rds"))
drcsmpls <- mipbigpanel$samples$Country == "DRC"

# Subset to DRC Samples
drcmips <- MIPanalyzer::filter_samples(x = mipbigpanel,
                                       sample_filter = drcsmpls,
                                       description = "Subset to DRC Samples")


# Correct for names in the adm1name in dhs recode
drcmips$samples$ADM1NAME <- gsub("Kongo Central", "Kongo-Central",
                                 drcmips$samples$ADM1NAME)

drcmips$samples$ADM1NAME <- gsub("Tanganyka", "Tanganyika",
                                           drcmips$samples$ADM1NAME)

saveRDS(drcmips, file = paste0(gdrive, "/data/derived_data/cd2013_dhs_bigbarcodemips.rds"))


# ghana focus for pulling apart the ibd
ghanasmpls <- mipbigpanel$samples$Country == "Ghana"

# Subset to DRC Samples
ghanamips <- MIPanalyzer::filter_samples(x = mipbigpanel,
                                       sample_filter = ghanasmpls,
                                       description = "Subset to Ghana Samples")

saveRDS(ghanamips, file = paste0(gdrive, "/data/derived_data/ghana_bigbarcodemips.rds"))

# ghana_drc focus for pulling apart the ibd
ghana_drcsmpls <- mipbigpanel$samples$Country == "Ghana" | mipbigpanel$samples$Country == "DRC"

# Subset to DRC Samples
ghana_drc_mips <- MIPanalyzer::filter_samples(x = mipbigpanel,
                                         sample_filter = ghana_drcsmpls,
                                         description = "Subset to Ghana and DRC Samples")


# Correct for names in the adm1name in dhs recode
ghana_drc_mips$samples$ADM1NAME <- gsub("Kongo Central", "Kongo-Central",
                                        ghana_drc_mips$samples$ADM1NAME)

ghana_drc_mips$samples$ADM1NAME <- gsub("Tanganyka", "Tanganyika",
                                        ghana_drc_mips$samples$ADM1NAME)


saveRDS(ghana_drc_mips, file = paste0(gdrive, "/data/derived_data/ghana_drc_bigbarcodemips.rds"))

