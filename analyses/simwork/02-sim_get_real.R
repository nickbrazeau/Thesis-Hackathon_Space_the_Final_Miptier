## .................................................................................
## Purpose: Batch these realized IBD calcs for drake
##
## Notes: Using "overly generous" IBD estimator from Verity, Aydemir, Brazeau et al.
##        End user may have to change number of nodes that they can parallelize across
##        for batch submission
## .................................................................................
library(tidyverse)
library(drake)
#............................................................
# read in and setup
#...........................................................
workers <- 150 # nodes to ask for, fewer nodes, less expensive for reading in data and not placing burden on scheduler
smpl_hosts <- readRDS("data/sim_data/sim_smpl_hosts.rds")

#............................................................
# functions for drake plan
#...........................................................
get_ibd_batch_wrap <- function(simdatpath, lvl, batchset, pair_hosts) {
  simdat <- readRDS(simdatpath)
  swfsim <- simdat %>%
    dplyr::filter(name == lvl)
  # pull out swfsim
  swfsim <- swfsim$swfsim[[1]]

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
# nest and split up pairwise comparisons for batching
#...........................................................
batchnum <- sort( rep(1:workers, ceiling(nrow(smpl_hosts) / workers)) )
batchnum <- batchnum[1 :nrow(smpl_hosts)]

smpl_hosts_nested <- smpl_hosts %>%
  dplyr::mutate(batchset = batchnum) %>%
  dplyr::group_by(batchset) %>%
  tidyr::nest() %>%
  dplyr::ungroup()

#............................................................
# now tidy up the rest in the inputs I need
#...........................................................
retmap <- expand.grid(c("mtn", "rift", "oppcorner"), 1:workers) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("lvl", "batchset")) %>%
  dplyr::mutate(simdatpath = "data/sim_data/swf_simulations.rds") %>%
  dplyr::left_join(., smpl_hosts_nested, by  = "batchset") %>%
  dplyr::rename(pair_hosts = data)

#............................................................
# Make Drake Plan
#...........................................................
batch_names <- paste0(retmap$lvl, retmap$batchset)
plan <- drake::drake_plan(
  sims = target(
    get_ibd_batch_wrap(simdatpath, lvl, batchset, pair_hosts),
    transform = map(
      .data = !!retmap,
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
     jobs = nrow(retmap),
     log_make = "swfsim_deploy_drake.log", verbose = 4,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     garbage_collection = TRUE,
     lock_envir = FALSE,
     lock_cache = FALSE)

