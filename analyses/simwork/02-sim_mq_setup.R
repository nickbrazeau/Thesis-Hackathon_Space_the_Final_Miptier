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
# QUESTION 1: Can we use Discent to do model comparisons/capture differences in geodistance
#...........................................................
#......................
# read in previous spatial setups
#......................
latticemodel <- readRDS("data/sim_data/lattice_model.rds")
gridmig <- readRDS("data/sim_data/gridmig_geodist.rds")
gridmigmat <- readRDS("data/sim_data/gridmig_prob_matrix.rds")


#............................................................
# run sWF simulator
#...........................................................
# magic numbers outside
Nesize <- 25
mscale <- 0.5
verity_coi2 <- readRDS("results/optim_coi_lambdas/optim_lambda.RDS")[2]

swf_sim_wrapper <- function(migmat, Nesize, mscale) {
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

  # from verity et al coi in the DRC: 2.23 (2.15– 2.31), going to use
  # a lambda of 2 as it is close
  # assuming deme size of 10 for ease
  # tlim at 10 generations as before from verity et al

  #......................
  # run structured WF
  #......................
  swfsim <- polySimIBD::sim_swf(pos =       pos,
                                migr_dist_mat = migmat,
                                N =         rep(Nesize, nrow(migmat)),
                                m =         rep(mscale, nrow(migmat)),
                                rho =       rho,
                                mean_coi =  rep(verity_coi2, nrow(migmat)),
                                tlim =      tlim)
  return(swfsim)
}

#......................
# run simulations
#......................
swfsim <- swf_sim_wrapper(gridmigmat, Nesize = Nesize, mscale = mscale)


#............................................................
# Pairwise IBD realizations
#...........................................................
# will assume 3 individuals in every deme
all_hosts <- tibble::tibble(host = 1:(Nesize*ncol(gridmigmat)),
                            host_deme = sort(rep(1:ncol(gridmigmat), Nesize)))
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
saveRDS(swfsim, "data/sim_data/mq_swf_simulations.rds")
saveRDS(smpl_hosts, "data/sim_data/sim_smpl_hosts_mq.rds")
