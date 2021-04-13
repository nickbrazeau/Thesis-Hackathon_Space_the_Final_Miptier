## .................................................................................
## Purpose: Simulation from the spatial DTDLsWF to gauge how
##          the discent model behaves under different "expected"
##          versus realized IBD pairings
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
library(polySimIBD)
set.seed(48)

#............................................................
# make spatial setup
#   plan will be a square always with some "genetic" barrier set of shapes
#   will use 'elevation' and euclidean distance btwn points to parameterize migration rates
#...........................................................
latticemodel <- readRDS("data/sim_data/lattice_model.rds") %>%
  dplyr::select(-c("migration", "deme")) # new migrations here; "new" deme names
simdat <- tibble::tibble(name = c("mtn", "rift", "oppcorner"),
                         gridmig = NA,
                         plotObj = NA)

#......................
# migration climbing central mountain
#......................
gridmig <- latticemodel %>%
  dplyr::mutate(migration = purrr::map2_dbl(longnum, latnum, function(x, y){
    mvtnorm::dmvnorm(c(x, y),
                     mean = c(50, 50),
                     sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                     log = T)

  }))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig))
raster::contour(raster::rasterFromXYZ(gridmig),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[1] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/250) )

#......................
# migration with central rift
#......................
gridmig <- latticemodel %>%
  dplyr::mutate(migration = purrr::map2_dbl(longnum, latnum, function(x, y){
    mvtnorm::dmvnorm(c(x, y),
                     mean = c(50, 50),
                     sigma = matrix(c(0.1, 1e-3, 1e-3, 5), ncol = 2),
                     log = T)

  }))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig))
raster::contour(raster::rasterFromXYZ(gridmig),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[2] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/250) )


#......................
# opposite corners
#......................
gridmig <- latticemodel %>%
  dplyr::mutate(migration = dplyr::case_when(
    longnum <= 50 & latnum <= 50 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      -mvtnorm::dmvnorm(c(x, y),
                        mean = c(25, 25),
                        sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                        log = T)}),

    longnum <= 50 & latnum > 50 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      mvtnorm::dmvnorm(c(x, y),
                       mean = c(25, 75),
                       sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                       log = T)}),
    longnum > 50 & latnum > 50 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      -mvtnorm::dmvnorm(c(x, y),
                        mean = c(75, 75),
                        sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                        log = T)}),

    longnum > 50 & latnum <= 50 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      mvtnorm::dmvnorm(c(x, y),
                       mean = c(75, 25),
                       sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                       log = T)})

  ))

# make other direction hill
gridmig_min <- min(gridmig$migration, na.rm = T)
gridmig_max <- max(gridmig$migration, na.rm = T)

gridmig <- gridmig %>%
  dplyr::mutate(migration = dplyr::case_when(
    longnum <= 50 & latnum > 50 ~ migration - gridmig_min,
    longnum <= 50 & latnum <= 50 ~ migration - gridmig_max,
    longnum > 50 & latnum > 50 ~ migration - gridmig_max,
    longnum > 50 & latnum <= 50 ~ migration - gridmig_min
  ))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig))
raster::contour(raster::rasterFromXYZ(gridmig),
                add = TRUE, drawlabels = FALSE, col = "#969696")





# store, same approx order of mag as distance for migration
simdat$gridmig[3] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/250) )




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
#   sample lattice locations
#...........................................................
locats <- latticemodel %>%
  dplyr::select(c("longnum", "latnum")) %>%
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
  locatcomb <- t(combn(locats$deme, 2)) %>%
    tibble::as_tibble(., .name_repair = "minimal") %>%
    magrittr::set_colnames(c("xy1", "xy2"))

  #......................
  # calculate
  #......................
  matlocat <- tibble::tibble(xy1 = locatcomb$xy1,
                             xy2 = locatcomb$xy2)

  matlocat$connval <- furrr::future_pmap_dbl(matlocat, function(xy1, xy2, locats){
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

  #......................
  # fix values in matrix
  #......................
  # spread out values for matrix
  matlocat_dist <- matlocat %>%
    tidyr::pivot_wider(data = .,
                       names_from = "xy2",
                       values_from = "connval")
  colnames(matlocat_dist)[1] <- locats$deme[1]
  matlocat_dist[,1] <- NA
  matlocat_dist <- rbind.data.frame(matlocat_dist, rep(NA, ncol(matlocat_dist)))
  rownames(matlocat_dist) <- locats$deme
  # convert to matrix
  matlocat_dist <- as.matrix(matlocat_dist)

  # make symmetrical
  matlocat_dist[lower.tri(matlocat_dist)]  <- t(matlocat_dist)[lower.tri(matlocat_dist)]
  diag(matlocat_dist) <- 0
  return(matlocat_dist)
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
    dplyr::mutate(migmat = purrr::map(distmat, function(x, scalar = 5){ # scaled above slightly too!
    x <- exp(-x/scalar)
    return(x)
  })
  )
# NB these are currently distances, later convert them into "connectivity"


# sanity
table( sample(1:ncol(simdat$migmat[[1]]), size = 1e3, prob = simdat$migmat[[1]][10,], replace = T) )
table( sample(1:ncol(simdat$migmat[[2]]), size = 1e3, prob = simdat$migmat[[2]][10,], replace = T) )
table( sample(1:ncol(simdat$migmat[[3]]), size = 1e3, prob = simdat$migmat[[3]][1,], replace = T) )



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
  # assuming deme size of 10 for ease
  # tlim at 10 generations as before from verity et al

  #......................
  # run structured WF
  #......................
  swfsim <- polySimIBD::sim_swf(pos =       pos,
                                migr_dist_mat = migmat,
                                N =         rep(25, nrow(migmat)),
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
all_hosts <- tibble::tibble(host = 1:(10*ncol(simdat$migmat[[1]])),
                            host_deme =sort(rep(1:ncol(simdat$migmat[[1]]), 10)))
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
dir.create("data/sim_data/")
# need to get correct "level" and unnest the list so we are swfsim class
saveRDS(simdat$swfsim[simdat$name == "mtn"][[1]], "data/sim_data/mtn_nonlinear_migration_sim_swf.rds")
saveRDS(simdat$swfsim[simdat$name == "rift"][[1]], "data/sim_data/rift_nonlinear_migration_sim_swf.rds")
saveRDS(simdat$swfsim[simdat$name == "oppcorner"][[1]], "data/sim_data/oppcorner_nonlinear_migration_sim_swf.rds")
# save out geodistances
saveRDS(simdat$distmat[simdat$name == "mtn"][[1]], "data/sim_data/mtn_nonlinear_migration_geodist.rds")
saveRDS(simdat$distmat[simdat$name == "rift"][[1]], "data/sim_data/rift_nonlinear_migration_geodist.rds")
saveRDS(simdat$distmat[simdat$name == "oppcorner"][[1]], "data/sim_data/oppcorner_nonlinear_migration_geodist.rds")
# save out sampled hosts
saveRDS(smpl_hosts, "data/sim_data/sim_smpl_hosts_nonlinear_migration.rds")
