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
# QUESTION 1: Can we get F values right based on deme/COI?
#...........................................................
#......................
# make spatial setup
#   plan will be a gradient of demes
#......................
nCell <- 100
coords <- round(seq(1, nCell, by = 10))
latticemodel <- expand.grid(coords, coords)
plot(latticemodel)
colnames(latticemodel) <- c("longnum", "latnum")
latticemodel <- latticemodel %>%
  dplyr::mutate(deme = 1:dplyr::n())

# euclidean for an iso by dist frame
# location combinatoins
locatcomb <- t(combn(sort(latticemodel$deme), 2)) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("deme1", "deme2"))

eucmigmat <- furrr::future_pmap_dbl(locatcomb, function(deme1, deme2){
  # get long lat
  xy1 <- latticemodel[latticemodel$deme == deme1, c("longnum", "latnum")]
  xy2 <- latticemodel[latticemodel$deme == deme2, c("longnum", "latnum")]

  # euclidean distance
  eucmigmat <- dist(rbind(xy1, xy2))
  return(as.numeric(eucmigmat))})

# save out for later
euc <- dplyr::bind_cols(locatcomb, distval = eucmigmat)
# spread out values for matrix
eucmigmat <- euc %>%
  tidyr::pivot_wider(data = .,
                     names_from = "deme2",
                     values_from = "distval")
eucmigmat[,1] <- NA
eucmigmat <- rbind.data.frame(eucmigmat, rep(NA, ncol(eucmigmat)))
# convert to matrix
eucmigmat <- as.matrix(eucmigmat)
# make symmetrical
eucmigmat[lower.tri(eucmigmat)]  <- t(eucmigmat)[lower.tri(eucmigmat)]
diag(eucmigmat) <- 0


#............................................................
# gradient
#...........................................................
coi_grad <- tibble::tibble(longnum = coords,
                           coigrad = seq(1, 8, length.out = length(coords))) %>%
  dplyr::left_join(latticemodel, ., by = "longnum") %>%
  dplyr::pull("coigrad")


ne_grad <- tibble::tibble(longnum = coords,
                          negrad = seq(5, 500, length.out = length(coords))) %>%
  dplyr::left_join(latticemodel, ., by = "longnum") %>%
  dplyr::pull("negrad")



#............................................................
# run sWF simulator
#...........................................................
swf_sim_wrapper <- function(migmat, coivec, nevec) {
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


  #......................
  # run structured WF
  #......................
  swfsim <- polySimIBD::sim_swf(pos =       pos,
                                migr_dist_mat = migmat,
                                N =         nevec,
                                m =         rep(0.5, nrow(migmat)),
                                rho =       rho,
                                mean_coi =  coivec,
                                tlim =      tlim)
  return(swfsim)
}

#............................................................
# run simulation for COI gradient
#...........................................................
coi_grad_sim <- swf_sim_wrapper(migmat = eucmigmat,
                                coivec = coi_grad,
                                nevec = rep(10, nrow(eucmigmat)))

# IBD realizations for coi grad
#   straightforward because host size doesn't change from 10
coi_all_hosts <- tibble::tibble(host = 1:(10*ncol(eucmigmat)),
                                host_deme = sort(rep(1:ncol(eucmigmat), 10)))
coi_all_hosts <- split(coi_all_hosts, factor(coi_all_hosts$host_deme))
coi_smpl_hosts <- lapply(coi_all_hosts, function(x){x[sample(1:nrow(x), size = 3),]}) %>%
  dplyr::bind_rows()

#......................
# expand out pairwise
#......................
coi_comb_hosts <- t(combn(1:nrow(coi_smpl_hosts), 2))
coi_comb_hosts <- split(coi_comb_hosts, 1:nrow(coi_comb_hosts))
exp_host_pairwise <- function(smpl_hosts, comb_hosts) {
  cbind.data.frame(smpl_hosts[comb_hosts[1], ],
                   smpl_hosts[comb_hosts[2], ]) %>%
    magrittr::set_colnames(c("smpl1", "deme1", "smpl2", "deme2"))
}
# get pairwise
coi_smpl_hosts <- lapply(coi_comb_hosts, exp_host_pairwise, smpl_hosts = coi_smpl_hosts) %>%
  dplyr::bind_rows()




#............................................................
# run simulation for Effective Population Size gradient
#............................................................
ne_grad_sim <- swf_sim_wrapper(migmat = eucmigmat,
                               coivec = rep(2.23, nrow(eucmigmat)),
                               nevec = ne_grad)

# will assume 3 individuals in every deme
ne_all_hosts <- tibble::tibble(host = 1:sum(ne_grad),
                               host_deme = rep(1:ncol(eucmigmat), time = ne_grad))
ne_all_hosts <- split(ne_all_hosts, factor(ne_all_hosts$host_deme))
ne_smpl_hosts <- lapply(ne_all_hosts, function(x){x[sample(1:nrow(x), size = 3),]}) %>%
  dplyr::bind_rows()

#......................
# expand out pairwise
#......................
ne_comb_hosts <- t(combn(1:nrow(ne_smpl_hosts), 2))
ne_comb_hosts <- split(ne_comb_hosts, 1:nrow(ne_comb_hosts))
exp_host_pairwise <- function(smpl_hosts, comb_hosts) {
  cbind.data.frame(smpl_hosts[comb_hosts[1], ],
                   smpl_hosts[comb_hosts[2], ]) %>%
    magrittr::set_colnames(c("smpl1", "deme1", "smpl2", "deme2"))
}
# get pairwise
ne_smpl_hosts <- lapply(ne_comb_hosts, exp_host_pairwise, smpl_hosts = ne_smpl_hosts) %>%
  dplyr::bind_rows()



#............................................................
# save out
#...........................................................
dir.create("data/sim_data/")
saveRDS(ne_grad_sim, "data/sim_data/ne_sim_grad.rds")
saveRDS(coi_grad_sim, "data/sim_data/coi_sim_grad.rds")
saveRDS(coi_smpl_hosts, "data/sim_data/sim_smpl_hosts_coi_fq.rds")
saveRDS(ne_smpl_hosts, "data/sim_data/sim_smpl_hosts_nepop_fq.rds")
saveRDS(euc, "data/sim_data/euclidean_geodist_fq.rds")
