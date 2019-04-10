#----------------------------------------------------------------------------------------------------
# Purpose of this script is to plot X-Y of space and genetic distance
#----------------------------------------------------------------------------------------------------
# libraries and imports
gdrive <- tcltk::tk_choose.dir()
library(tidyverse)

drcmips <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_gen_space_epi_final.rds"))


summary(drcmips$results$admin_gen_dist$IBS_verity_maskhets$adim_gen_dist)
summary(drcmips$results$admin_gen_dist$IBS_verity_majallele$adim_gen_dist)
summary(drcmips$results$admin_gen_dist$DAB_malariagen_wsaf$adim_gen_dist)
summary(drcmips$results$admin_gc_dist$distance)

dplyr::left_join(x = drcmips$results$admin_gen_dist$IBS_verity_majallele,
                 y = drcmips$results$admin_gc_dist,
                 by = c("admin1_1", "admin1_2")) %>%
  ggplot() +
  geom_point(aes(x= distance, y = adim_gen_dist)) +
  scale_x_log10() +
  theme_bw()



