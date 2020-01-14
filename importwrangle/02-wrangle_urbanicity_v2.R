#..............................................................
# Purpose of this script is to wrangle urbanicity in and
# define cities. Going to use kmeans clustering
#..............................................................
library(raster)
library(tidyverse)
source("R/basics.R")

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
# Need to Scale Data for Kmeans
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
# Perform kmeans
#..............................................................
https://rspatial.org/raster/rs/4-unsupclassification.html
rasters <- stack(c(trav, fric, nightlights, worldpop))
rasters.pts <- values(rasters)
rasters.pts.full <- na.omit(rasters.pts)

#.............
# K-means Clustering of DRC for urbanicity
#.............
keda <- data.frame(k = seq(2, 200, by = 1))
keda$kmeans <- map(keda$k,
                   function(k){
                     return(kmeans(x =rasters.pts.full, centers = k))})
keda$wss <- map(keda$kmeans, "withinss")
keda$totalwss <- map(keda$wss, function(x){return(sum(x))})


keda.df <- keda %>%
  dplyr::select(c("k", "totalwss")) %>%
  tidyr::unnest(cols = totalwss)

kesteda.plotObj <- keda.df %>%
  tibble::as.tibble(.) %>%
  ggplot() +
  geom_line(aes(x=k, y=totalwss)) +
  geom_point(aes(x=k, y=totalwss)) +
  geom_vline(xintercept = 15, color = "red", linetype = 2, alpha = 0.8) +
  theme_minimal() +
  ylab("Total Within-Cluster Sum of Squares") +
  xlab("K")



#.............
# PCA on random sample
#.............
pca <- stats::prcomp(rasters.pts)

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
urbanicity <- raster::predict(rasters,
                              pca,
                              index=1) # just use PC1

# invert scale, PCA has cities going towards negative
values(urbanicity) <- values(urbanicity) * -1



#................
# Out
#................
dir.create("data/derived_data/urbanicity_raster")
raster::writeRaster(urbanicity,
                    filename = "data/derived_data/urbanicity_raster/urbanicity.grd")


#..............................................................
# Urabnicity as a Catchment Area
#..............................................................
# Urbanitity in the top third quartile is "city-ish" and below is rural.
urban.class <- urban
thirdquart <- quantile(values(urban), probs = 0.75, na.rm = T)
values(urban.class)[values(urban.class) < thirdquart] <- 0
values(urban.class)[values(urban.class) >= thirdquart] <- 1
ggplot() +
  ggspatial::layer_spatial(data = urban, aes(fill = stat(band1))) +
  scale_fill_distiller("Urbanicity Score", type = "div", palette = "RdYlBu") +
  prettybasemap_nodrc +
  theme(legend.position = "right",
        legend.title = element_text(size = 13, face = "bold", angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 12, face = "bold", angle = 0))

raster::writeRaster(urban.class,
                    filename = "data/derived_data/urbanicity_raster/urbanicity_binary.grd")



