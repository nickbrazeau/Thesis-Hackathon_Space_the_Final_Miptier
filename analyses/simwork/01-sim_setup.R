## .................................................................................
## Purpose: Simulation from the spatial DTDLsWF to gauge "power" of discent model
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
#remotes::install_github("nickbrazeau/polySimIBD")
library(polySimIBD)

#............................................................
# make spatial setup
#   plan will be a square always with some "genetic" barrier set of shapes
#   will use 'elevation' and euclidean distance btwn points to parameterize migration rates
#...........................................................
set.seed(48)
nCell <- 300
simdat <- tibble::tibble(name = c("mtn", "rift", "fourcorner"),
                         gridmig = NA,
                         plotObj = NA)

#......................
# migration climbing central mountain
#......................
gridmig <- matrix(NA, nrow = nCell, ncol = nCell)
gridmig <- data.frame(longnum = as.vector(row(gridmig)),
                      latnum = as.vector(col(gridmig)),
                      migration = NA)
gridmig <- gridmig %>%
  dplyr::mutate(migration = purrr::map2_dbl(longnum, latnum, function(x, y){
    mvtnorm::dmvnorm(c(x, y),
                     mean = c(150, 150),
                     sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                     log = T)

  }))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig))
raster::contour(raster::rasterFromXYZ(gridmig),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[1] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/1e3) )

#......................
# migration with central rift
#......................
gridmig <- matrix(NA, nrow = nCell, ncol = nCell)
gridmig <- data.frame(longnum = as.vector(row(gridmig)),
                      latnum = as.vector(col(gridmig)),
                      migration = NA)
gridmig <- gridmig %>%
  dplyr::mutate(migration = purrr::map2_dbl(longnum, latnum, function(x, y){
    mvtnorm::dmvnorm(c(x, y),
                     mean = c(150, 150),
                     sigma = matrix(c(0.1, 1e-3, 1e-3, 5), ncol = 2),
                     log = T)

  }))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig))
raster::contour(raster::rasterFromXYZ(gridmig),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[2] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/1e3) )


#......................
# four corners
#......................
gridmig <- matrix(NA, nrow = nCell, ncol = nCell)
gridmig <- data.frame(longnum = as.vector(row(gridmig)),
                      latnum = as.vector(col(gridmig)),
                      migration = NA)
gridmig <- gridmig %>%
  dplyr::mutate(migration = dplyr::case_when(
    longnum <= 100 & latnum <= 100 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      mvtnorm::dmvnorm(c(x, y),
                       mean = c(50, 50),
                       sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                       log = T)}),
    longnum <= 100 & latnum > 200 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      mvtnorm::dmvnorm(c(x, y),
                       mean = c(50, 250),
                       sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                       log = T)}),

    longnum <= 300 & latnum <= 100  &  longnum > 200  ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      mvtnorm::dmvnorm(c(x, y),
                       mean = c(250, 50),
                       sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                       log = T)}),
    longnum <= 300  &  longnum > 200 & latnum > 200 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      mvtnorm::dmvnorm(c(x, y),
                       mean = c(250, 250),
                       sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                       log = T)}),

  ))

# make "ground"
gridmig$migration[is.na(gridmig$migration)] <- quantile(gridmig$migration,
                                                        probs = 0.01,
                                                        na.rm = T)
# visualize to confirm
plot(raster::rasterFromXYZ(gridmig))
raster::contour(raster::rasterFromXYZ(gridmig),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[3] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/1e2) )




#......................
# store nice plots
#......................
simdat <- simdat %>%
  dplyr::mutate(plotObj = purrr::map(gridmig, function(x){
    x %>%
      ggplot() +
      geom_tile(aes(x = longnum, y = latnum, fill = migration)) +
      geom_contour(aes(x = longnum, y = latnum, z = migration), color = "black") +
      scale_fill_viridis_c("Migration Topology") +
      coord_fixed() +
      theme_void()
  }))

#............................................................
# Location Sampling and migration rate calculations
#   sample 350 locations
#...........................................................
locats <- gridmig[sample(1:nrow(gridmig), size = 350), c("longnum", "latnum")] %>%
  dplyr::mutate(deme = 1:dplyr::n())

get_dist_matrix  <- function(gridmig, locats) {
  #......................
  # liftovers
  #......................
  # storage mat
  mat <- matrix(NA, nrow = nrow(locats), ncol = nrow(locats))
  # convert to raster
  rstr <- raster::rasterFromXYZ(gridmig)
  #......................
  # get combinations I need
  #......................
  locatcomb <- t(combn(1:nrow(locats), 2)) %>%
    tibble::as_tibble(., .name_repair = "minimal") %>%
    magrittr::set_colnames(c("xy1", "xy2"))

  #......................
  # calculate
  #......................
  matlocat <- tibble::tibble(xy1 = locatcomb$xy1,
                             xy2 = locatcomb$xy2)

  connval <- furrr::future_pmap(matlocat, function(xy1, xy2, locats){
    # get long lat
    xy1 <- locats[xy1, c("longnum", "latnum")]
    xy2 <- locats[xy2, c("longnum", "latnum")]

    # get height
    height1 <- raster::extract(rstr, xy1)
    height2 <- raster::extract(rstr, xy2)

    # euclidean distance
    euc <- dist(rbind(xy1, xy2))
    # "connectedness" is 1/distance plus difference in "heights" -- i.e. if same "plane" or not
    height <- sqrt((height2 - height1)^2)
    ret <- euc +  height
    return(ret)},
    locats = locats)

  # spread out values for matrix
  for (i in 1:nrow(locatcomb)) {
    x <- as.numeric(attr(connval[[i]], "Labels")[1])
    y <- as.numeric(attr(connval[[i]], "Labels")[2])
    mat[x,y] <- connval[[i]]
  }

  # make symmetrical
  mat[lower.tri(mat)]  <- t(mat)[lower.tri(mat)]
  diag(mat) <- 0
  return(mat)
}

#......................
# find distance matrices
#......................
simdat <- simdat %>%
  dplyr::mutate(distmat = purrr::map(gridmig, get_dist_matrix, locats = locats))
# NB these are currently distances, later convert them into "connectivity"

#......................
# liftover to migration matrix
#......................
simdat <- simdat %>%
  dplyr::mutate(migmat = purrr::map(distmat, function(x, scalar = 1e3){
                                      x <- exp(-x/scalar)
                                      return(x)
                                    })
  )
# NB these are currently distances, later convert them into "connectivity"




#............................................................
# run sWF simulator
#...........................................................
swf_sim_wrapper <- function(migmat) {
  #......................
  # magic numbers
  #......................
  # Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
  # Aimee gets this number by taking the inverse of Mile's estimate of the CO recombination rate of 13.5 kb/cM
  rho <- 7.4e-7
  # going to assume we can only detect things 10 generations ago
  tlim <- 10

  # approximate average of Pf3d7 Chromosome Lengths
  pflen <- 1.664e6
  # assuming single chromosome for ease
  # assuming 1e3 loci
  pos <- sort(sample(1.664e6, 1e3))

  # from verity et al coi in the DRC: 2.23 (2.15– 2.31)
  # assuming deme size of 50 for ease
  # tlim at 10 generations as before from verity et al

  #......................
  # run structured WF
  #......................
  swfsim <- polySimIBD::sim_swf(pos =       pos,
                                migr_dist_mat = migmat,
                                N =         rep(50, nrow(migmat)),
                                m =         rep(0.5, nrow(migmat)),
                                rho =       rho,
                                mean_coi =  rep(2.23, nrow(migmat)),
                                tlim =      10)
  return(swfsim)
}

#......................
# run simulations
#......................
simdat <- simdat %>%
  dplyr::mutate(swfsim = purrr::map(migmat, swf_sim_wrapper))


#............................................................
# Pairwise IBD realizations
#...........................................................
# will assume 3 individuals in every deme
all_hosts <- tibble::tibble(host = 1:(50*ncol(simdat$migmat[[1]])),
                            host_deme =sort(rep(1:ncol(simdat$migmat[[1]]), 50)))
all_hosts <- split(all_hosts, factor(all_hosts$host_deme))
smpl_hosts <- lapply(all_hosts, function(x){x[sample(1:nrow(x), size = 3),]}) %>%
  dplyr::bind_rows()

#......................
# expand out pairwise
#......................
comb_hosts <- t(combn(1:nrow(smpl_hosts), 2))
comb_hosts <- split(comb_hosts, 1:nrow(comb_hosts))
exp_host_pairwise <- function(smpl_hosts, comb_hosts) {
  cbind.data.frame(smpl_hosts[comb_hosts[1], ],
                   smpl_hosts[comb_hosts[2], ]) %>%
    magrittr::set_colnames(c("smpl1", "deme1", "smpl2", "deme2"))
}
# get pairwise
smpl_hosts <- lapply(comb_hosts, exp_host_pairwise, smpl_hosts = smpl_hosts) %>%
  dplyr::bind_rows()

#......................
# save out
#......................
saveRDS(simdat, "data/sim_data/swf_simulations.rds")
locats %>%
  dplyr::mutate(deme = 1:dplyr::n()) %>%
  saveRDS(., "data/sim_data/sim_locations.rds")
saveRDS(smpl_hosts, "data/sim_data/sim_smpl_hosts.rds")
