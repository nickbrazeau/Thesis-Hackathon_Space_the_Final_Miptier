#................................................................................................................................
# Purpose of this script is to recode urbanicity to better
# reflect DRC boundaries/continuum
#................................................................................................................................

# IMPORTS and dependencies
library(tidyverse)
library(sf)
library(plotly)
source("R/basics.R")
tol <- 1e-3
set.seed(48)

# Notes on Lonely PSUs
# http://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
options(survey.lonely.psu="adjust")

# spatial from the DHS -- these are cluster level vars
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")  %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))
sf::st_geometry(dt) <- NULL
#.............
# weights
#.............
dt <- dt %>%
  dplyr::mutate(hv005_wi = hv005/1e6
  )


ge <- sf::st_as_sf(readRDS(file = "data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
ge <- ge %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))


#..............................
#### A note on urbanicity ####
#..............................

# Potential (significant) misclassification bias in the DHS DRC-II coding of
# urban vs. rural as has been noted here https://journals.sagepub.com/doi/10.1177/0021909617698842
# and can be seen by comparing hv025/026 with population density, light density, build, etc.

# #.............
# # Urban from DHS
# #.............
# levels(factor(haven::as_factor(dt$hv025))) # no missing
# dt$hv025_fctb <- haven::as_factor(dt$hv025)
# dt$hv025_fctb = forcats::fct_relevel(dt$hv025_fctb, "rural")
# # check
# xtabs(~hv025 + hv025_fctb, data = dt, addNA = T)


#.............
# cluster degree of "build"
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2014
summary(dt$built_population_2014)
hist(dt$built_population_2014)
sum(dt$built_population_2014 < 0.01)
median( dt$built_population_2014[dt$built_population_2014 < 0.05] )
hist( dt$built_population_2014[dt$built_population_2014 < 0.05] )
# DECISION: Will use a logit transformation to get back to the real-line (and scale)
# large number of 0s
dt$built_population_2014_scale <- my.scale(logit(dt$built_population_2014, tol = tol), center = T, scale = T) # use logit to transform to real line
hist(dt$built_population_2014_scale)
summary(dt$built_population_2014_scale); sd(dt$built_population_2014_scale) # despite skew, scale seems to work

#.............
# cluster night-time light density
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(dt$nightlights_composite)
hist(dt$nightlights_composite)
hist( dt$nightlights_composite[dt$nightlights_composite < 0.05] )
hist( dt$nightlights_composite[dt$nightlights_composite > 0.05] )
# large number of 0s (again)
dt$nightlights_composite_scale <- my.scale(log(dt$nightlights_composite + tol), center = T, scale = T)
hist(dt$nightlights_composite_scale)
summary(dt$nightlights_composite_scale); sd(dt$nightlights_composite_scale) # despite skew, scale seems to work


#.............
# cluster worldpop population-density estimate
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(dt$all_population_count_2015)
hist(dt$all_population_count_2015)
hist( dt$all_population_count_2015[dt$all_population_count_2015 < 5e4] )
hist( dt$all_population_count_2015[dt$all_population_count_2015 > 5e4] )
# no 0s here but a lot of small pops
dt$all_population_count_2015_scale <- my.scale(log(dt$all_population_count_2015 + tol), center = T, scale = T)
hist(dt$all_population_count_2015_scale)
summary(dt$all_population_count_2015_scale); sd(dt$all_population_count_2015_scale) # scale here seems to compensate


#.............
# Accessibility from Weiss
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(dt$travel_times_2015)
hist(dt$travel_times_2015)
dt <- dt %>%
  dplyr::mutate(travel_times_2015_scale = my.scale(log(travel_times_2015 + tol), center = T, scale = T))

hist(dt$travel_times_2015_scale) # many, many 0s -- these are urban centers/places near big towns
hist(dt$travel_times_2015_scale) # standardization doesn't look as good, but should capture urban v. rural well
summary(dt$travel_times_2015_scale); sd(dt$travel_times_2015_scale)




#.............
# PCA
#.............
# compute PCA, drop to single observations (since clusters all have some obs)
# Otherwise inflate our degree of certainty in our PCA
urbanmat <- dt[!duplicated(dt$hv001),
               c("hv001", "hv025", "built_population_2014_scale", "nightlights_composite_scale",
                 "all_population_count_2015_scale", "travel_times_2015_scale")]

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
pcaeigplot <- tibble(rururb = haven::as_factor(urbanmat$hv025),
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
# K-means Clustering of Urbanicity
#.............
# Apply k-means with k=4
k <- kmeans(urbanmatpca$x[,1:4], 2, nstart=25, iter.max=1000)
k # the explanation of between SS versus total SS is ~65.9%


# plot new designations
pcaeig_newrural <- tibble(rururb_new = factor(k$cluster),
                          rururb_orig = haven::as_factor(urbanmat$hv025),
                          pc1 = urbanmatpca$x[,1], pc2 =urbanmatpca$x[,2], pc3 = urbanmatpca$x[,3])

# seems reasonable for what I want, a lot of the originally coded urban centers
# are recoded to rural, which is more in line with what we would expect
ggplot() +
  geom_jitter(data = pcaeig_newrural, aes(x=pc1, y = pc2, shape = rururb_orig, color = rururb_new))



#.............
# Degree of Urbanicity
#.............
urbanicity <- cbind(urbanmat[,c("hv001", "hv025")],
                    urbanpcascore_cont_scale_clst = urbanmatpca$x[,1],
                    urban_rural_pca_fctb_clst = factor(k$cluster, levels = c(1,2), labels = c("rural", "urban"))
)
urbanicity$urban_rural_pca_fctb_clst <- relevel(urbanicity$urban_rural_pca_fctb_clst, "urban")


urbanicity.plot <- urbanicity %>%
  left_join(x = ., y = ge, by = "hv001")

ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(data = urbanicity.plot,
             aes(x = longnum, y = latnum,
                 color = urbanpcascore_cont_scale_clst,
                 shape = urban_rural_pca_fctb_clst)) +
  viridis::scale_color_viridis("viridis")


#.............
# write out
#.............
saveRDS(urbanicity, file = "data/derived_data/DHS_urbanicity_recode.rds")
