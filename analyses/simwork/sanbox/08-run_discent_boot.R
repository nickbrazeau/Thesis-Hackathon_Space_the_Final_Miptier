## .................................................................................
## Purpose: Run discent many times to get a "bootstrapped" CI of the gradient descent F and M estimates
## Notes:
## .................................................................................
#............................................................
#  setup
#...........................................................
library(tidyverse)
library(drake)

#......................
# drank function
#......................
discent_wrapper <- function(rep, q, f_start, m_start, learn) {
  qbang <- enquo(q)
  input <- readRDS("data/sim_data/sim_gengeodat.rds") %>%
    dplyr::filter(q == !!qbang)

  input <- input$gengeodat[[1]] %>%
    dplyr::filter(locat1 != locat2)
  our_start_params <- rep(f_start, 100) # 100 sim clusters
  names(our_start_params) <- 1:100
  our_start_params <- c(our_start_params, "m" = m_start)
  ret <- discent::deme_inbreeding_spcoef(K_gendist_geodist = input,
                                         start_params = our_start_params,
                                         m_lowerbound = -.Machine$double.xmax,
                                         m_upperbound = 100,
                                         learningrate = learn,
                                         steps = 1e4,
                                         report_progress = FALSE,
                                         return_verbose = FALSE)
  # now write out
  dir.create("/pine/scr/n/f/nfb/Projects/Space_the_Final_Miptier/finalboots_sim_cluster_inbreed_ests/", recursive = T)
  saveRDS(ret,
          file = paste0("/pine/scr/n/f/nfb/Projects/Space_the_Final_Miptier/finalboots_sim_cluster_inbreed_ests/",
                        "bootsimdat_", q, "_rep", rep,  ".RDS"))
  return(0)
}

#......................
# read in best results
#......................
beststart <- readRDS("results/simulated_min_cost_inbreedingresults/sim_min_cost_inbreedingresults.RDS")
beststart <- beststart %>%
  dplyr::select(c("q", "f_start", "m_start", "learn"))

# expand out this param tables
reps <- 100
repbeststart <- lapply(1:reps, function(x){return(beststart)}) %>%
  dplyr::bind_rows(.id = "rep") %>%
  dplyr::arrange(q, rep)

#......................
# make drake plan
#......................
rep_names <- paste0("rep", repbeststart$q, repbeststart$rep)
plan <- drake::drake_plan(
  runs = target(
    discent_wrapper(rep, q, f_start, m_start, learn),
    transform = map(
      .data = !!repbeststart,
      .names = !!rep_names
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
     jobs = nrow(repbeststart),
     log_make = "finalboot_simwork_discent_drc_dat_deploy_drake.log", verbose = 4,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     garbage_collection = TRUE,
     lock_envir = FALSE,
     lock_cache = FALSE)
