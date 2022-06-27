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
# QUESTION 2: Can we get F values right based on deme/COI?
#...........................................................
#......................
# read in spatial setup and then make gradients: (1) COI grad; (2) Ne grad
#   plan will be a gradient of demes (using previous euclidean setup)
#......................
eucmigmat <- readRDS("data/sim_data/euclidean_prob_matrix.rds")
latticemodel <- readRDS("data/sim_data/lattice_model.rds")
coords <- round(seq(1, nrow(latticemodel), by = 10)) # just need moving along x-axis


#............................................................
# gradient
#...........................................................
optim_lambda_from_verity_etal <- readRDS("results/optim_coi_lambdas/optim_lambda.RDS")
# coi
coi_grad <- tibble::tibble(longnum = coords,
                           coigrad = seq(min(optim_lambda_from_verity_etal),
                                         max(optim_lambda_from_verity_etal),
                                         length.out = length(coords))) %>%
  dplyr::left_join(latticemodel, ., by = "longnum")
# sanity
plot(raster::rasterFromXYZ(xyz = coi_grad[c("longnum", "latnum", "coigrad")]))
# pull
coi_grad <- coi_grad %>% dplyr::pull("coigrad")

# ne
ne_grad <- tibble::tibble(longnum = coords,
                          negrad =  round(seq(5, 50, length.out = length(coords)))) %>%
  dplyr::left_join(latticemodel, ., by = "longnum")
# sanity
plot(raster::rasterFromXYZ(xyz = ne_grad[c("longnum", "latnum", "negrad")]))
# pull
ne_grad <- ne_grad %>%
  dplyr::pull("negrad")



#............................................................
# run sWF simulator
#...........................................................
# magic numbers outside
Nesize <- 25
mscale <- 0.5

swf_sim_wrapper <- function(migmat, coivec, nevec, mscale) {
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
                                m =         rep(mscale, nrow(migmat)),
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
                                nevec = rep(Nesize, nrow(eucmigmat)),
                                mscale = mscale)

# IBD realizations for coi grad
#   straightforward because host size doesn't change from Nesize
coi_all_hosts <- tibble::tibble(host = 1:(Nesize*ncol(eucmigmat)),
                                host_deme = sort(rep(1:ncol(eucmigmat), Nesize)))
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
                               coivec = rep(optim_lambda_from_verity_etal[2], nrow(eucmigmat)),
                               nevec = ne_grad,
                               mscale = mscale)

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
# saveRDS(euc, "data/sim_data/euclidean_geodist_fq.rds") remember, euclidean is same as in mq setup
