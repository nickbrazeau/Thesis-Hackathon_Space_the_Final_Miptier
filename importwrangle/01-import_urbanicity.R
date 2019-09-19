#................................................................................................................................
# Purpose of this script is to recode urbanicity based on several
# collinear cluster variables. Collinearity can be handled by some
# methods but not by others (e.g. glmms). So we will use a PCA
# to extract out the signal
#................................................................................................................................

library(tidyverse)
library(sf)
library(plotly)
tol <- 1e-3
source("R/basics.R")



DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))
ge <- sf::st_as_sf(readRDS(file = "data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
ge <- ge %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

gc <- readRDS("data/raw_data/dhsdata/datasets/CDGC62FL.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(hv001 = dhsclust)

ge_gc <- left_join(ge, gc, by = "hv001")

#..............................
# A note on urbanicity
#..............................
# Going to use four variables to get at urbanicity: Build, night-light density, worldpop densities, and travel times

#.............
# cluster degree of "build"
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2014
summary(ge_gc$built_population_2014)
hist(ge_gc$built_population_2014)
sum(ge_gc$built_population_2014 < 0.01)
median( ge_gc$built_population_2014[ge_gc$built_population_2014 < 0.05] )
hist( ge_gc$built_population_2014[ge_gc$built_population_2014 < 0.05] )
# DECISION: Will use a logit transformation to get back to the real-line (and scale)
# large number of 0s
ge_gc$built_population_2014_scale <- my.scale(logit(ge_gc$built_population_2014, tol = tol), center = T, scale = T) # use logit to transform to real line
hist(ge_gc$built_population_2014_scale)
summary(ge_gc$built_population_2014_scale); sd(ge_gc$built_population_2014_scale) # despite skew, scale seems to work

#.............
# cluster night-time light density
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(ge_gc$nightlights_composite)
hist(ge_gc$nightlights_composite)
hist( ge_gc$nightlights_composite[ge_gc$nightlights_composite < 0.05] )
hist( ge_gc$nightlights_composite[ge_gc$nightlights_composite > 0.05] )
# large number of 0s (again)
ge_gc$nightlights_composite_scale <- my.scale(log(ge_gc$nightlights_composite + tol), center = T, scale = T)
hist(ge_gc$nightlights_composite_scale)
summary(ge_gc$nightlights_composite_scale); sd(ge_gc$nightlights_composite_scale) # despite skew, scale seems to work


#.............
# cluster worldpop population-density estimate
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(ge_gc$all_population_count_2015)
hist(ge_gc$all_population_count_2015)
hist( ge_gc$all_population_count_2015[ge_gc$all_population_count_2015 < 5e4] )
hist( ge_gc$all_population_count_2015[ge_gc$all_population_count_2015 > 5e4] )
# no 0s here but a lot of small pops
ge_gc$all_population_count_2015_scale <- my.scale(log(ge_gc$all_population_count_2015 + tol), center = T, scale = T)
hist(ge_gc$all_population_count_2015_scale)
summary(ge_gc$all_population_count_2015_scale); sd(ge_gc$all_population_count_2015_scale) # scale here seems to compensate


#.............
# Accessibility from Weiss
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(ge_gc$travel_times_2015)
hist(ge_gc$travel_times_2015)
ge_gc <- ge_gc %>%
  dplyr::mutate(travel_times_2015_scale = my.scale(log(travel_times_2015 + tol), center = T, scale = T))

hist(ge_gc$travel_times_2015_scale) # many, many 0s -- these are urban centers/places near big towns
hist(ge_gc$travel_times_2015_scale) # standardization doesn't look as good, but should capture urban v. rural well
summary(ge_gc$travel_times_2015_scale); sd(ge_gc$travel_times_2015_scale)




#.............
# PCA
#.............
# compute PCA, drop to single observations (since clusters all have some obs)
# Otherwise inflate our degree of certainty in our PCA
urbanmat <- ge_gc[,
               c("hv001", "urban_rura", "built_population_2014_scale", "nightlights_composite_scale",
                 "all_population_count_2015_scale", "travel_times_2015_scale")]
sf::st_geometry(urbanmat) <- NULL

urbanmatpca <- prcomp(urbanmat[, c("built_population_2014_scale", "nightlights_composite_scale",
                                   "all_population_count_2015_scale", "travel_times_2015_scale")])

# compute variance explained and Zi values
urbanmatpca$var <- (urbanmatpca$sdev ^ 2) / sum(urbanmatpca$sdev ^ 2) * 100
urbanmatpca$loadings <- abs(urbanmatpca$rotation)
urbanmatpca$loadings <- sweep(urbanmatpca$loadings, 2, colSums(urbanmatpca$loadings), "/")


#.............
# Analyze PCA Seperation
#.............
# make a plot
pcaeigplot <- tibble(rururb = haven::as_factor(urbanmat$urban_rura),
                     pc1 = urbanmatpca$x[,1], pc2 =urbanmatpca$x[,2], pc3 = urbanmatpca$x[,3])

plotly::plot_ly(pcaeigplot, x = ~pc1, y = ~pc3, z = ~pc2, color = ~rururb)%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC3'),
                      zaxis = list(title = 'PC2')))
#.............
# View on Biplot
#.............
library(ggfortify)
ggplot2::autoplot(urbanmatpca,
                  loadings = TRUE,
                  loadings.label = TRUE,
                  loadings.label.size  = 3)


#.............
# Degree of Urbanicity
#.............
ge_gc <- ge_gc %>%
  dplyr::mutate(urbanpcascore = urbanmatpca$x[,1])


ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(data = ge_gc,
             aes(x = longnum, y = latnum,
                 color = urbanpcascore,
                 shape = haven::as_factor(urban_rura))) +
  viridis::scale_color_viridis("viridis")




#.............
# write out
#.............
ret <- ge_gc %>%
  dplyr::select(c("hv001", "built_population_2014_scale", "nightlights_composite_scale",
                  "all_population_count_2015_scale", "travel_times_2015_scale", "urbanpcascore"))

sf::st_geometry(ret) <- NULL

saveRDS(ret, file = "data/derived_data/DRC_urbanicity.rds")
