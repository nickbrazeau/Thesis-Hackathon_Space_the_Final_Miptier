## .................................................................................
## Purpose: Wrangle Mobility data
##
## Notes:
## .................................................................................
library(mobility)
library(tidyverse)
library(raster)
#............................................................
# Wrangle Population Data Data
#...........................................................
# worldpop already cropped to DRC from worldpop site
worldpopraw <- raster::raster("data/raw_data/worldpop/cod_ppp_2013.tif")
# increase size for tract
worldpopraw <- raster::aggregate(worldpopraw, fact = 10, fun = sum)

#......................
# voroni territories popN
#......................
vrdf <- readRDS("data/distance_data/voroni_base.RDS")
voroni_popN <- rep(NA, nrow(vrdf))
for (i in 1:nrow(vrdf)){
  voroni_popN[i] <- unlist( raster::extract(x = worldpopraw,
                                            y = sf::as_Spatial(vrdf$x[i]),
                                            fun = sum,
                                            na.rm = T) ) # unlist for matrix
}

#......................
# voroni territories popN
#......................
vrdf <- readRDS("data/distance_data/voroni_base.RDS")
voroni_popN <- rep(NA, nrow(vrdf))
for (i in 1:nrow(vrdf)){
  voroni_popN[i] <- unlist( raster::extract(x = worldpopraw,
                                            y = sf::as_Spatial(vrdf$x[i]),
                                            fun = sum,
                                            na.rm = T) ) # unlist for matrix
}

# names
names(voroni_popN) <- vrdf$IPUMSID

#......................
# provinces popN
#......................
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds")
prov_popN <- rep(NA, nrow(DRCprov))
for (i in 1:nrow(DRCprov)){
  prov_popN[i] <- unlist( raster::extract(x = worldpopraw,
                                          y = sf::as_Spatial(DRCprov$geometry[i]),
                                          fun = sum,
                                          na.rm = T) ) # unlist for matrix
}
# add names
names(prov_popN) <- DRCprov$adm1name
prov_popN <- round(prov_popN)

#............................................................
# Build Mobility Matrices
#...........................................................
#......................
# voroni
#......................
voroni_flows <- readRDS("data/distance_data/voroni_migration_flows_fromworldpop.RDS")

# Build mobility matrix from the longform data
M_voroni <- get_mob_matrix(orig = voroni_flows$NODEI,
                           dest = voroni_flows$NODEJ,
                           value = voroni_flows$PrdMIG)

#......................
# anq123
#......................
anq_flows <- readRDS("data/distance_data/anq123_migration_flows.RDS")

# Build mobility matrix from the longform data
M_anq <- get_mob_matrix(orig = anq_flows$org,
                        dest = anq_flows$dest,
                        value = anq_flows$mgrtn)


#............................................................
# Build Distance Matrices
#...........................................................
#......................
# voroni
#......................
vrdf_longlat <- tibble::tibble(IPUMSID = vrdf$IPUMSID)
vrdf_longlat$geometry <- sf::st_centroid(vrdf$x)
vrdf_longlat <- sf::st_sf(vrdf_longlat)
sf::st_crs(vrdf_longlat) <- 4326 # assume planar

D_voroni <- sf::st_distance(vrdf_longlat, which = "Great Circle")
# remove units https://stackoverflow.com/questions/46935207/removing-units-from-an-r-vector
attr(D_voroni, "units") <- NULL
class(D_voroni) <- setdiff(class(D_voroni),"units")
colnames(D_voroni) <- vrdf_longlat$IPUMSID
rownames(D_voroni) <- vrdf_longlat$IPUMSID


#......................
# anq123
#......................
anq_longlat <- tibble::tibble(adm1name = DRCprov$adm1name)
anq_longlat$geometry <- sf::st_centroid(DRCprov$geometry)
anq_longlat <- sf::st_sf(anq_longlat)
sf::st_crs(anq_longlat) <- 4326 # assume planar

D_anq <- sf::st_distance(anq_longlat, which = "Great Circle")
# remove units for easier downstream
# https://stackoverflow.com/questions/46935207/removing-units-from-an-r-vector
attr(D_anq, "units") <- NULL
class(D_anq) <- setdiff(class(D_anq),"units")
colnames(D_anq) <- anq_longlat$adm1name
rownames(D_anq) <- anq_longlat$adm1name

#............................................................
# Fitting Mobility Model
#...........................................................
#......................
# voroni
#......................
# make assumption of diagonal being zero for fit
diag(M_voroni) <- 0
voroni_total <- rowSums(M_voroni, na.rm=TRUE)
prob_trav <- summarize_mobility(
  fit_prob_travel(travel = voroni_total - diag(M_voroni), total = voroni_total)
)


#......................
# anq123
#......................
# make assumption of diagonal being zero for fit
anq_total <- rowSums(M_anq, na.rm=TRUE)
prob_trav <- summarize_mobility(
  fit_prob_travel(travel = anq_total - diag(M_anq), total = anq_total)
)
mod_anq <- summarize_mobility(
  fit_mobility(M_anq, D_anq, prov_popN), ac_lags = 5
)




#............................................................
# check goodnes of fit
#...........................................................
#......................
# anq123
#......................
check_mobility(M = M_anq,
               D = D_anq,
               N = prov_popN,
               mod = mod_anq)

mod_sim_anq <- sim_mobility(D = D_anq,
                            N = prov_popN,
                            mod = mod_anq,
                            n = 3)
mod_sim_anq[1:5,1:5]



ggplot(data=melt(mod_sim_anq)) +
  geom_tile(aes(x=factor(Var2),
                y=factor(Var2),
                fill=value)) +
  xlab('Destination') + ylab("Origin") +
  theme_bw() + theme(axis.text.x=element_text(size=10),
                     axis.text.y=element_text(size=10),
                     axis.title.x=element_text(size=12, margin = margin(t = 15)),
                     axis.title.y=element_text(size=12, margin = margin(r = 15)),
                     legend.position='bottom') +
  viridis::scale_fill_viridis(option='inferno', direction=1) +
  guides(fill=guide_colorbar(title='Observed trips',
                             title.position='top',
                             label.theme=element_text(size=9),
                             barwidth=20,
                             barheight=0.5,
                             frame.colour='black',
                             ticks=TRUE))
