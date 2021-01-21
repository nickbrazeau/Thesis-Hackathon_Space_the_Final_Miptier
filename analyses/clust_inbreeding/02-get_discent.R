## .................................................................................
## Purpose: Setup the the Malecot's Spatial Gradient Descent (discent) model and run it with drake on cluster
##
## Notes:
## .................................................................................
#............................................................
#  setup
#...........................................................
workers <- 512 # nodes to ask for, fewer nodes, less expensive for reading in data and not placing burden on scheduler
library(tidyverse)
library(drake)

# cluster names
mtdt <- readRDS("data/derived_data/sample_metadata.rds")
clsts <- sort(unique(c(mtdt$hv001)))


#............................................................
# functions for drake
#...........................................................
drake_wrapper <- function(batchset_df) {

  # call future
  no_cores <- future::availableCores() - 1
  if (no_cores > 1) {
    future::plan(future::multicore, workers = no_cores)
  } else {
    future::plan("sequential")
  }

  #......................
  # internal function to wrap discent
  #......................
  discent_wrapper <- function(inputpath, f_start, m_start, f_learn, m_learn, clst_names, id) {
    input <- readRDS(as.character(inputpath))
    our_start_params <- rep(f_start, 351) # 351 is all clusters in drc
    names(our_start_params) <- clst_names
    our_start_params <- c(our_start_params, "m" = m_start)
    ret <- discent::deme_inbreeding_spcoef(K_gendist_geodist = input,
                                           start_params = our_start_params,
                                           m_lowerbound = 1e-25,
                                           m_upperbound = 5e-4,
                                           f_learningrate = f_learn,
                                           m_learningrate = m_learn,
                                           full_matrix = FALSE,
                                           steps = 1e4,
                                           report_progress = FALSE,
                                           return_verbose = FALSE)
    return(ret)
  }

  #......................
  # run batches to not overload scheduler
  #......................
  batchset_df$discentret = furrr::future_pmap(batchset_df, discent_wrapper)

  # now write out
  dir.create("/pine/scr/n/f/nfb/Projects/Space_the_Final_Miptier/cluster_inbreed_ests/", recursive = T)
  saveRDS(batchset_df,
          file = paste0("/pine/scr/n/f/nfb/Projects/Space_the_Final_Miptier/cluster_inbreed_ests/",
                        "gendat_", unique(batchset_df$id), ".RDS")
  )
  return(0)
}

#............................................................
# make the parammap we are going to explore
#...........................................................
allsmpls_gengeodatpaths <- list.files("data/derived_data/allsmpls_clst_inbreeding_dat/", full.names = T)
monoclonals_gengeodatpaths <- list.files("data/derived_data/coione_clst_inbreeding_dat/", full.names = T)
gengeodatpaths <- c(allsmpls_gengeodatpaths, monoclonals_gengeodatpaths)
fs <- seq(0.1, 0.9, by = 0.2)
ms <- c(1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5)
f_learningrate <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)
m_learningrate <- c(1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12)

param_map <- expand.grid(gengeodatpaths, fs, ms, f_learningrate, m_learningrate) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("inputpath", "f_start", "m_start", "f_learn", "m_learn")) %>%
  dplyr::mutate(clst_names = list(clsts))


#............................................................
# make the batches and drake plan
#...........................................................
#......................
# nest and split up pairwise comparisons for batching
#......................
batchnum <- sort( rep(1:workers, ceiling(nrow(param_map) / workers)) )
batchnum <- batchnum[1 :nrow(param_map)]

param_map_nested <- param_map %>%
  dplyr::mutate(id = batchnum,
                batchset = batchnum) %>%
  dplyr::group_by(batchset) %>%
  tidyr::nest() %>%
  dplyr::ungroup()

#......................
# make drake plan
#......................
batch_names <- paste0("batch", param_map_nested$batchset)
plan <- drake::drake_plan(
  runs = target(
    drake_wrapper(data),
    transform = map(
      .data = !!param_map_nested,
      .names = !!batch_names
    )
  )
)

#......................
# call drake to send out to slurm
#......................
options(clustermq.scheduler = "slurm",
        clustermq.template = "drake_clst/slurm_clustermq_LL_long_litemem.tmpl")
make(plan,
     parallelism = "clustermq",
     jobs = nrow(param_map_nested),
     log_make = "discent_drc_dat_deploy_drake.log", verbose = 4,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     garbage_collection = TRUE,
     lock_envir = FALSE,
     lock_cache = FALSE)


