## .................................................................................
## Purpose: Determine computational time of discent on normal macbook
##
## Notes:
## .................................................................................
#............................................................
#  setup
#...........................................................
library(tidyverse)
library(drake)
library(discent)

#............................................................
# Simulation function
#...........................................................
#' @description not a generalizable function -- but to make a quick wrapper
#' @param nsmpls integer; samples per deme
#' @param ndemes integer; number of demes

time_elapsed_discent <- function(nsmpls = 10, ndemes = 10, steps = 100, rep = 1) {
  #......................
  # framework setup and tidy
  #......................
  locats <- as.character(sample(1:.Machine$integer.max, size = ndemes, replace = F))
  locats <- sort(rep(locats, times = nsmpls))
  locats <- tibble::tibble(smpl = paste0("smpl", as.character(sample(1:.Machine$integer.max, size = length(locats), replace = F))),
                           locat = locats)
  # bring together
  K_gendist_geodist <- as.data.frame(t(combn(locats$smpl, 2)))
  colnames(K_gendist_geodist)[1] <- "smpl1"
  colnames(locats)[1] <- "smpl1"
  K_gendist_geodist <- dplyr::left_join(K_gendist_geodist, locats, by = "smpl1")
  colnames(K_gendist_geodist)[2] <- "smpl2"
  colnames(locats)[1] <- "smpl2"
  K_gendist_geodist <- dplyr::left_join(K_gendist_geodist, locats, by = "smpl2")
  colnames(K_gendist_geodist) <- c("smpl1", "smpl2", "locat1", "locat2")

  # now remove selfs (not necessary since we sample w/out replacement above, but for transparency)
  K_gendist_geodist <- K_gendist_geodist %>%
    dplyr::filter(locat1 != locat2)

  #......................
  # simulate genetic and geographic distance
  #......................
  # make up some genetic distance
  K_gendist_geodist$gendist <- rbeta(nrow(K_gendist_geodist), 0.5, 0.5)

  geodist <- K_gendist_geodist %>%
    dplyr::select(c("locat1", "locat2")) %>%
    dplyr::filter(!duplicated(.))

  # make up some geographic distance
  geodist <- geodist %>%
    dplyr::mutate(geodist = rlnorm(dplyr::n(), 3, 1))

  #......................
  # tidy and start params
  #......................
  # locat names
  locatnames <- unique(c(geodist$locat1, geodist$locat2))
  # bring together
  K_gendist_geodist <- dplyr::left_join(K_gendist_geodist, geodist, by = c("locat1", "locat2"))

  start_params = rep(0.5, length(locatnames))
  names(start_params) <- locatnames
  start_params <- c(start_params, "m" = 1e-2)
  #......................
  # run model
  #......................
  start <- Sys.time()
  mod <- discent::deme_inbreeding_spcoef(K_gendist_geodist = K_gendist_geodist,
                                         start_params = start_params,
                                         m_lowerbound = 0,
                                         m_upperbound = 1,
                                         learningrate = 1e-10,
                                         steps = steps,
                                         report_progress = FALSE)
  time_elapsed <- Sys.time() - start
  # save out
  dat <- tibble::tible(
    nsmpls = nsmpls,
    ndemes = ndemes,
    steps = steps,
    time_elapsed = time_elapsed)
  # now write out
  dir.create("results/discent_comp_time/", recursive = T)
  saveRDS(batchset_df,
          file = paste0("results/discent_comp_time/",
                        "comptime_", paste0(nsmpls, "-", ndemes, "-", steps, "-", rep), ".RDS")
  )
  return(0)

}


#............................................................
# Param Map for Testing
#...........................................................
nsmpls <- c(1, 3, 5, 10, 25, 50)
ndemes <- round(seq(1, 1000, length.out = 25))
steps <- c(100, 1e3, 1e4, 1e5)
param_map <- expand.grid(nsmpls, ndemes, steps) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("nsmpls", "ndemes", "steps"))
# make reps
nreps <- 5
param_map <- lapply(1:nreps, function(x)return(param_map))
param_map <- param_map %>%
  dplyr::bind_rows(., .id = "rep")

#......................
# make drake plan
#......................
run_names <- paste0("run", 1:nrow(param_map))
plan <- drake::drake_plan(
  runs = target(
    time_elapsed_discent(nsmpls, ndemes, steps, rep),
    transform = map(
      .data = !!param_map,
      .names = !!run_names
    )
  )
)



make(plan,
     log_make = "computational_time_discent.log", verbose = 4,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     garbage_collection = TRUE,
     lock_envir = FALSE,
     lock_cache = FALSE)
