## .................................................................................
## Purpose: Setup the the Malecot's Spatial Gradient Descent (discent) model and run it with drake on cluster
##
## Notes:
## .................................................................................
#............................................................
#  setup
#...........................................................
workers <- 8000 # slurm array jobs to partition across
library(tidyverse)
library(drake)


#............................................................
# functions for drake
#...........................................................
drake_wrapper <- function(batchset_df, batchset) {

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
  discent_wrapper <- function(inputpath, f_start, m_start, learn) {
    input <- readRDS(as.character(inputpath)) %>%
      dplyr::filter(locat1 != locat2)
    # cluster details
    #   cluster are different sizes depending on coi of 1 or all
    clstnum <- length(unique(c(input$locat1, input$locat2)))
    clst_names <- sort(unique(c(input$locat1, input$locat2)))
    # start param
    our_start_params <- rep(f_start, clstnum)
    names(our_start_params) <- clst_names
    our_start_params <- c(our_start_params, "m" = m_start)
    ret <- discent::deme_inbreeding_spcoef(K_gendist_geodist = input,
                                           start_params = our_start_params,
                                           m_lowerbound = -.Machine$double.xmax,
                                           m_upperbound = 100,
                                           f_learningrate = learn,
                                           m_learningrate = learn,
                                           steps = 5e4,
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
                        "gendat_", unique(batchset), ".RDS")
  )
  return(0)
}

#............................................................
# make the parammap we are going to explore
#...........................................................
allsmpls_gengeodatpaths <- list.files("data/derived_data/allsmpls_clst_inbreeding_dat/", full.names = T)
monoclonals_gengeodatpaths <- list.files("data/derived_data/coione_clst_inbreeding_dat/", full.names = T)
gengeodatpaths <- c(allsmpls_gengeodatpaths, monoclonals_gengeodatpaths)
fs <- seq(0, 1, by = 0.1)
ms <- c(1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)
learningrate <- c(1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,
                  1e-9, 1e-8, 1e-7, 1e-6, 1e-5)
param_map <- expand.grid(gengeodatpaths, fs, ms, learningrate) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("inputpath", "f_start", "m_start", "learn"))


#............................................................
# make the batches and drake plan
#...........................................................
#......................
# nest and split up pairwise comparisons for batching
#......................
batchnum <- sort( rep(1:workers, ceiling(nrow(param_map) / workers)) )
batchnum <- batchnum[1 :nrow(param_map)]

param_map_nested <- param_map %>%
  dplyr::mutate(batchset = batchnum) %>%
  dplyr::group_by(batchset) %>%
  tidyr::nest() %>%
  dplyr::ungroup()

#......................
# make drake plan
#......................
batch_names <- paste0("batch", param_map_nested$batchset)
plan <- drake::drake_plan(
  runs = target(
    drake_wrapper(data, batchset),
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


