## .................................................................................
## Purpose: Simulation from the spatial DTDLsWF to gauge "power" of discent model
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
#remotes::install_github("nickbrazeau/polySimIBD", ref = "develop")
library(polySimIBD)
set.seed(48)


#............................................................
# QUESTION 1: Can we compare model distance matrices 
#...........................................................
#......................
# make spatial setup
#   plan will be a square always with some "mountain" barrier
#......................
nCell <- 100
coords <- round(seq(1, nCell, by = 10))
latticemodel <- expand.grid(coords, coords)
plot(latticemodel)
colnames(latticemodel) <- c("longnum", "latnum")


#......................
# migration climbing central mountain
#......................
latticemodel <- latticemodel %>%
  dplyr::mutate(migration = purrr::map2_dbl(longnum, latnum, function(x, y){
    mvtnorm::dmvnorm(c(x, y),
                     mean = c(nCell/2, nCell/2),
                     sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
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
diag(gridmigmat) <- 0


#......................
# calculate for euclidean
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
diag(eucmat) <- 0


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
                                N =         rep(10, nrow(migmat)),
                                m =         rep(0.5, nrow(migmat)),
                                rho =       rho,
                                mean_coi =  rep(2.23, nrow(migmat)),
                                tlim =      tlim)
  return(swfsim)
}

#......................
# run simulations
#......................
swfsim <- swf_sim_wrapper(gridmigmat)


#............................................................
# Pairwise IBD realizations
#...........................................................
# will assume 3 individuals in every deme
all_hosts <- tibble::tibble(host = 1:(10*ncol(gridmigmat)),
                            host_deme = sort(rep(1:ncol(gridmigmat), 10)))
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
saveRDS(swfsim, "data/sim_data/mq_swf_simulations.rds")
saveRDS(gridmig, "data/sim_data/gridmig_geodist.rds")
saveRDS(euc, "data/sim_data/euclidean_geodist.rds")
saveRDS(smpl_hosts, "data/sim_data/sim_smpl_hosts_mq.rds")
