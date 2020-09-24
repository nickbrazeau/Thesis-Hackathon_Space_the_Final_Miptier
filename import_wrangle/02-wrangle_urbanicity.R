#..............................................................
# Purpose of this script is to wrangle urbanicity in to a single
# measure for analysis. This needs to be done to account for the
# fact that night light density is pretty sparse in the DRC (only very major cities)
# as well as to account for potential extreme levels of collinearity in for population
# density, travel times, etc.
#..............................................................
library(raster)
library(tidyverse)
library(sf)
source("R/basics.R")
set.seed(48)
# create bounding box of Central Africa for space
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"
#..............................................................
# Read in Data
#..............................................................
travraw <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_accessibility_to_cities_v1.0_latest_10_.18_40_8_2020_09_18.tiff")
fricraw <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_friction_surface_v1_Decompressed_latest_10_.18_40_8_2020_09_18.tiff")
travraw <- raster::crop(travraw, caf)
fricraw <- raster::crop(fricraw, caf)

# nightlights cropped from earlier
nightlightsraw <- raster::raster("data/derived_data/nightlights/drc_nightlights_raw.grd")
# worldpop already cropped to DRC from worldpop site
worldpopraw <- raster::raster("data/raw_data/worldpop/cod_ppp_2013.tif")


# sanity
raster::crs(travraw)
raster::crs(fricraw)
raster::crs(nightlightsraw)
raster::crs(worldpopraw)

#......................
# misc
#......................
# take care of missing values
values(travraw)[values(travraw) == -9999] <- NA
# need to apply masking early on for nightlights due to extreme values (for not land, and corrections, etc)
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
raster::crs(DRC)
nightlightsraw <- raster::mask(nightlightsraw, mask = DRC)

#..............................................................
# Take to Same Spatial Scale
#..............................................................
# agg up
nightlightsraw <- raster::aggregate(nightlightsraw, fact = 2, fun = mean)
worldpopraw <- raster::aggregate(worldpopraw, fact = 10, fun = mean)

#..........................................................
# Let's make same resolution and extent
# will use travel as common factor
#..........................................................
urbnrstr <- list(travraw, fricraw, nightlightsraw, worldpopraw)
urbnrstr.res <- lapply(urbnrstr, function(x){
  x <- raster::projectRaster(x, travraw) # to make all same extent
  return(x)
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
values(fric)[!is.na(values(fric))] <- my.scale(logit(values(fric)[!is.na(values(fric))], tol = 0.1))
sd(values(fric), na.rm = T)

# nightlights
summary(values(nightlights))
hist(values(nightlights))
# these less than zero values arose from tranformation, drop
values(nightlights)[values(nightlights) < 0 ] <- NA

# transformation does not result because of extreme number of zeros. Variacnce should be ok
nightlights <- raster::scale(nightlights, center = T, scale = T)
sd(values(nightlights), na.rm = T)

# world pop
summary(values(worldpop))
hist(values(worldpop))

# again a lot of zeroes
worldpop <- raster::scale(worldpop, center = T, scale = T)
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
summary(values(urbanicity))

#................
# reproject
#................
urbanicity.reproj <- raster::projectRaster(urbanicity, crs = raster::crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs"))


#................
# Out
#................
dir.create("data/derived_data/urbanicity_raster")
# reproj
raster::writeRaster(urbanicity.reproj,
                    filename = "data/derived_data/urbanicity_raster/urbanicity_reproj.grd",
                    overwrite = T)

# basic
raster::writeRaster(urbanicity,
                    filename = "data/derived_data/urbanicity_raster/urbanicity.grd",
                    overwrite = T)



# Save out raw
raster::writeRaster(urbanicity.raw,
                    filename = "data/derived_data/urbanicity_raster/urbanicity_raw_nottrunc.grd",
                    overwrite = T)

