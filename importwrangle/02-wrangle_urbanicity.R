#..............................................................
# Purpose of this script is to wrangle urbanicity in to a single
# measure for analysis. This needs to be done to account for the
# fact that night light density is pretty spare in the DRC (only very major cities)
# as well as to account for potential extreme levels of collinearity in the
# data for the spatial random forest
#..............................................................
library(raster)
library(tidyverse)
source("R/basics.R")
set.seed(48)
# need this for bounding box
DRC <- sf::as_Spatial(osmdata::getbb("Democratic Republic of the Congo",
                                     featuretype = "country",
                                     format_out = 'sf_polygon'))

#..............................................................
# Read in Data
#..............................................................
trav <- raster::raster("data/derived_data/MAPrasters/trav_acc.gri")
fric <- raster::raster("data/derived_data/MAPrasters/frict.grd")
nightlights <- raster::raster("data/derived_data/nightlights/drc_nightlights_raw.grd")
worldpop <- raster::raster("data/derived_data/worldpop/drc_worldpop_raw.grd")

#..............................................................
# Take to Same Spatial Scale
#..............................................................
# agg up
trav <- raster::aggregate(trav, fact = 6, fun = mean)
fric <- raster::aggregate(fric, fact = 6, fun = mean)

nightlights <- raster::aggregate(nightlights, fact = 12, fun = mean)
worldpop <- raster::aggregate(worldpop, fact = 60, fun = mean)

#..........................................................
# Let's make same resolution and extent
# will use travel as common factor
#..........................................................
urbnrstr <- list(trav, fric, nightlights, worldpop)
urbnrstr.res <- lapply(urbnrstr, function(x){

  ret <- raster::mask(x, DRC)

  ret <- raster::projectRaster(x,
                               trav)
  return(ret)
})

trav <- urbnrstr.res[[1]]
fric <- urbnrstr.res[[2]]
nightlights <- urbnrstr.res[[3]]
worldpop <- urbnrstr.res[[4]]


#..............................................................
# Need to Scale Data for PCA
#..............................................................
# travel
summary(values(trav))
hist(values(trav))
# log distributed
values(trav) <- my.scale(log(values(trav) + 0.1))
sd(values(trav), na.rm = T)

# fric
summary(values(fric))
hist(values(fric))
# log distributed
# transformation does not result in normal but will have to do
values(fric) <- my.scale(logit(values(fric), tol = 0.1))
sd(values(fric), na.rm = T)

# nightlights
summary(values(nightlights))
hist(values(nightlights))
# these less than zero values arose from tranformation, drop
values(nightlights)[values(nightlights) < 0 ] <- NA

# transformation does not result because of extreme number of zeros. Variacnce should be ok
nightlights <-raster::scale(nightlights, center = T, scale = T)
sd(values(nightlights), na.rm = T)

# world pop
summary(values(worldpop))
hist(values(worldpop))
# these less than zero values arose from tranformation, drop
values(worldpop)[values(worldpop) < 0 ] <- NA

# again a lot of zeroes
worldpop <-raster::scale(worldpop, center = T, scale = T)
sd(values(worldpop), na.rm = T)


#..............................................................
# Perform PCA
#..............................................................
# https://stackoverflow.com/questions/19866009/pca-using-raster-datasets-in-r
rasters <- stack(c(trav, fric, nightlights, worldpop))
names(rasters) <- c("trav", "fric", "nightlight", "worldpop")
rasters.pts <- values(rasters)
rasters.pts.boolean <- apply(rasters.pts, 1, function(x) return(!any(is.na(x))))

#.............
# PCA on random sample
#.............
pca <- stats::prcomp(rasters.pts[rasters.pts.boolean,])

# compute variance explained and Zi values
pca$var <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2) * 100
pca$loadings <- abs(pca$rotation)
pca$loadings <- sweep(pca$loadings, 2, colSums(pca$loadings), "/")

# look at var explained
pca$var


#.............
# Analyze PCA Seperation
#.............
library(plotly)
pcaeigplot <- tibble(pc1 = pca$x[,1],
                     pc2 =pca$x[,2],
                     pc3 = pca$x[,3])

plotly::plot_ly(pcaeigplot, x = ~pc1, y = ~pc3, z = ~pc2)%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC3'),
                      zaxis = list(title = 'PC2')))
#.............
# View on Biplot
#.............
library(ggfortify)
ggplot2::autoplot(pca,
                  loadings = TRUE,
                  loadings.label = TRUE,
                  loadings.label.size  = 3)

#..............................................................
# Predict New Raster Surface from PCA
#..............................................................
urbanicity.raw <- urbanicity <- raster::predict(rasters,
                                                pca,
                                                index=1) # just use PC1
# cap it at the 99.9% for outliers
val <- quantile(values(urbanicity), c(0.999), na.rm = T)
values(urbanicity)[values(urbanicity) >= val] <- val
plot(urbanicity)



#................
# Out
#................
dir.create("data/derived_data/urbanicity_raster")
raster::writeRaster(urbanicity,
                    filename = "data/derived_data/urbanicity_raster/urbanicity.grd",
                    overwrite = T)

#..............................................................
# Save out raw
#..............................................................
raster::writeRaster(urbanicity.raw,
                    filename = "data/derived_data/urbanicity_raster/urbanicity_raw_nottrunc.grd",
                    overwrite = T)

#..............................................................
# Urbanicity as a Binary
#..............................................................
urbanicity.binary <- raster::predict(rasters,
                                     pca,
                                     index=1) # just use PC1
urbanicity.binary.bool <- is.na(values(urbanicity.binary))
cutoffval <- quantile(urbanicity.binary, 0.95, na.rm = T)

values(urbanicity.binary)[!urbanicity.binary.bool] <-
  as.numeric( values(urbanicity.binary)[!urbanicity.binary.bool] > cutoffval )

plot(urbanicity.binary)

#..............................................................
# out
#..............................................................
raster::writeRaster(urbanicity.binary,
                    filename = "data/derived_data/urbanicity_raster/urbanicity_binary.grd",
                    overwrite = T)

