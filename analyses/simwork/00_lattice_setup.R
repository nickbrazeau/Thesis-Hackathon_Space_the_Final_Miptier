## .................................................................................
## Purpose: Make a lattice model to use for the polySimIBD modeling framework
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
set.seed(48)

#............................................................
# Make a spatial setup to ask two question (01; 02 scripts)
#     Will use a Lattice Model
#...........................................................
#......................
# Lattice
#......................
nCell <- 100
coords <- round(seq(1, nCell, by = 10))
latticemodel <- expand.grid(coords, coords)
plot(latticemodel)
colnames(latticemodel) <- c("longnum", "latnum")


#......................
# migration with mountain peak
#......................
latticemodel <- latticemodel %>%
  dplyr::mutate(migration = purrr::map2_dbl(longnum, latnum, function(x, y){
    mvtnorm::dmvnorm(c(x, y),
                     mean = c(nCell/2, nCell/2),
                     sigma = matrix(c(0.1, 0, 0, 0.1), ncol = 2),
                     log = T)
  }))

# visualize to confirm
plot(raster::rasterFromXYZ(latticemodel))
raster::contour(raster::rasterFromXYZ(latticemodel),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
latticemodel <- latticemodel %>%
  dplyr::mutate(migration = migration/1e2,
                deme = 1:dplyr::n())


#............................................................
# true distance matrix based on gridmig
#...........................................................
# convert to raster
rstr <- raster::rasterFromXYZ(latticemodel)

#......................
# get combinations I need
#......................
locatcomb <- t(combn(sort(values(rstr$deme)), 2)) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("deme1", "deme2"))

#......................
# calculate for euclidean and height
#......................
# euc and height
gridmig <- furrr::future_pmap_dbl(locatcomb, function(deme1, deme2){
  # get long lat
  xy1 <- latticemodel[latticemodel$deme == deme1, c("longnum", "latnum")]
  xy2 <- latticemodel[latticemodel$deme == deme2, c("longnum", "latnum")]

  # get height
  height1 <- raster::extract(rstr, xy1)[, "migration"]
  height2 <- raster::extract(rstr, xy2)[, "migration"]

  # euclidean distance
  euc <- dist(rbind(xy1, xy2))
  # "connectedness" is 1/distance plus difference in "heights" -- i.e. if same "plane" or not
  height <- sqrt((height2 - height1)^2)
  ret <- euc +  height
  return(ret)})

# save out for later
gridmig <- dplyr::bind_cols(locatcomb, distval = gridmig)
# spread out values for matrix
gridmigmat <- gridmig %>%
  tidyr::pivot_wider(data = .,
                     names_from = "deme2",
                     values_from = "distval")
gridmigmat[,1] <- NA
gridmigmat <- rbind.data.frame(gridmigmat, rep(NA, ncol(gridmigmat)))
# convert to matrix
gridmigmat <- as.matrix(gridmigmat)
# make symmetrical
gridmigmat[lower.tri(gridmigmat)]  <- t(gridmigmat)[lower.tri(gridmigmat)]

# liftover to migration PROBABILITY matrix w/out the diagonal
gridmigmat <- exp(-gridmigmat)
gridmigmat <- gridmigmat/rowSums(gridmigmat, na.rm = T)

# now add in diagonal and offset for stay more often -- 50% time at home
max_val <- max(gridmigmat, na.rm = T)
diag(gridmigmat) <- 1
# now re-make probability matrix
gridmigmat <- gridmigmat/rowSums(gridmigmat)

# sanity
table( sample(1:ncol(gridmigmat), size = 1e3, prob = gridmigmat[1,], replace = T) )
table( sample(1:ncol(gridmigmat), size = 1e3, prob = gridmigmat[56,], replace = T) )

#......................
# calculate for euclidean
#   NB -- this euclidean matrix is not use in this script -- it is saved for use in the Fq script
#......................
# euclidean
euc <- furrr::future_pmap_dbl(locatcomb, function(deme1, deme2){
  # get long lat
  xy1 <- latticemodel[latticemodel$deme == deme1, c("longnum", "latnum")]
  xy2 <- latticemodel[latticemodel$deme == deme2, c("longnum", "latnum")]

  # euclidean distance
  euc <- dist(rbind(xy1, xy2))
  return(as.numeric(euc))})

# save out for later
euc <- dplyr::bind_cols(locatcomb, distval = euc)
# spread out values for matrix
eucmat <- euc %>%
  tidyr::pivot_wider(data = .,
                     names_from = "deme2",
                     values_from = "distval")
eucmat[,1] <- NA
eucmat <- rbind.data.frame(eucmat, rep(NA, ncol(eucmat)))
# convert to matrix
eucmat <- as.matrix(eucmat)
# make symmetrical
eucmat[lower.tri(eucmat)]  <- t(eucmat)[lower.tri(eucmat)]

# liftover to migration PROBABILITY matrix w/out the diagonal
eucmat <- exp(-eucmat)
eucmat <- eucmat/rowSums(eucmat, na.rm = T)

# now add in diagonal element
max_val <- max(eucmat, na.rm = T)
# some offset for more stay -- about 50% of time want them to stay
diag(eucmat) <- 1

# liftover to migration PROBABILITY matrix WITH the diagonal
eucmat <- eucmat/rowSums(eucmat)

# sanity
table( sample(1:ncol(eucmat), size = 1e3, prob = eucmat[1,], replace = T) )
table( sample(1:ncol(eucmat), size = 1e3, prob = eucmat[56,], replace = T) )

# quick compare
plot(euc$distval, gridmig$distval)

#............................................................
# Save out
#...........................................................
dir.create("data/sim_data/")
saveRDS(latticemodel, "data/sim_data/lattice_model.rds")
saveRDS(eucmat, "data/sim_data/euclidean_prob_matrix.rds")
saveRDS(gridmigmat, "data/sim_data/gridmig_prob_matrix.rds")
saveRDS(gridmig, "data/sim_data/gridmig_geodist.rds")
saveRDS(euc, "data/sim_data/euclidean_geodist.rds")
