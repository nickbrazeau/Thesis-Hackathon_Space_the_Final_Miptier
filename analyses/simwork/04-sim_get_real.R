## .................................................................................
## Purpose: Batch these realized IBD calcs for drake
##
## Notes: Using "overly generous" IBD estimator from Verity, Aydemir, Brazeau et al.
##        End user may have to change number of nodes that they can parallelize across
##        for batch submission
## .................................................................................
library(tidyverse)
library(drake)
workers <- 8000
#............................................................
# functions for drake plan
#...........................................................
get_ibd_batch_wrap <- function(simdatpath, lvl, batchset, pair_hosts) {
  # readin swfsim
  swfsim <- readRDS(simdatpath)

  # unnest pairwairse comparisons
  pair_hosts$simIBD <- pmap_dbl(pair_hosts[,c("smpl1", "smpl2")], get_ver_realized_ibd, swfsim = swfsim)

  # now write out
  dir.create("/pine/scr/n/f/nfb/Projects/Space_the_Final_Miptier/swf_sim_ret/", recursive = T)
  saveRDS(pair_hosts,
          file = paste0("/pine/scr/n/f/nfb/Projects/Space_the_Final_Miptier/swf_sim_ret/",
                        lvl, "_", "batchset_", batchset, ".RDS"
                        )
          )
  return(0)

}

#......................
# get verity realized ibd
#......................
get_ver_realized_ibd <- function(swfsim, smpl1, smpl2) {
  # pull out ARG pieces
  # trees point left --  don't care about number of strains w/in, if any ibd, then loci is ibd
  arg <- polySimIBD::get_arg(swf = swfsim, host_index = c(smpl1, smpl2))
  conn <- purrr::map(arg, "c")

  # now look at loci by loci
  get_loci_pairwise_ibd <- function(conni, this_coi) {
    smpl1con <- conni[1:this_coi[1]]
    smpl2con <- conni[(this_coi[1]+1):(cumsum(this_coi)[2])]
    # get IBD connections between 1 and 2
    pwconn <- which(smpl2con %in% 0:(this_coi[1]-1) )
    # if any, IBD
    pwconn <- sum(pwconn)
    return(pwconn >= 1)
  }
  # iterate through
  lociibd <- purrr::map_dbl(.x = conn,
                             .f = get_loci_pairwise_ibd, this_coi = swfsim$coi[c(smpl1, smpl2)])
  # loci IBD
  return(sum(lociibd)/length(lociibd))
}


#............................................................
# read in data and make retmap
#...........................................................
mq_smpl_hosts <- readRDS("data/sim_data/sim_smpl_hosts_mq.rds")
coi_smpl_hosts <- readRDS("data/sim_data/sim_smpl_hosts_coi_fq.rds")
ne_smpl_hosts <- readRDS("data/sim_data/sim_smpl_hosts_nepop_fq.rds")

# bring together
retmap <- tibble::tibble(
  lvl = c("mq", "coi", "ne"),
  simdatpath = c("data/sim_data/mq_swf_simulations.rds",
                 "data/sim_data/coi_sim_grad.rds",
                 "data/sim_data/ne_sim_grad.rds"),
  hosts = list(mq_smpl_hosts,
               coi_smpl_hosts,
               ne_smpl_hosts))


#............................................................
# nest and split up pairwise comparisons for batching
#...........................................................
batchnum <- sort( rep(1:workers, ceiling( sum(sapply(retmap$hosts, nrow)) / workers )) )
batchnum <- batchnum[1 : sum(sapply(retmap$hosts, nrow))]

retmap_nested <- retmap %>%
  tidyr::unnest(cols = "hosts") %>%
  dplyr::mutate(batchset = batchnum) %>%
  dplyr::group_by(batchset, lvl, simdatpath) %>%
  tidyr::nest() %>%
  dplyr::ungroup() %>%
  dplyr::rename(pair_hosts = data)

#............................................................
# Make Drake Plan
#...........................................................
batch_names <- paste0(retmap_nested$lvl, retmap_nested$batchset)
plan <- drake::drake_plan(
  sims = target(
    get_ibd_batch_wrap(simdatpath, lvl, batchset, pair_hosts),
    transform = map(
      .data = !!retmap_nested,
      .names = !!batch_names
    )
  )
)


#......................
# call drake to send out to slurm
#......................
options(clustermq.scheduler = "slurm",
        clustermq.template = "drake_clst/slurm_clustermq_LL.tmpl")
make(plan,
     parallelism = "clustermq",
     jobs = nrow(retmap_nested),
     log_make = "swfsim_deploy_drake.log", verbose = 4,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     garbage_collection = TRUE,
     lock_envir = FALSE,
     lock_cache = FALSE)

