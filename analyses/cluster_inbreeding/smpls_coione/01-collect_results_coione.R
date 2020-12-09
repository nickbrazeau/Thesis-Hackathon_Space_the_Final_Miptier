#########################################################################
# Purpose: Collect the results from the Malecot's Spatial Gradient Descent (cluster level)
#
# Notes:
#########################################################################
library(tidyverse)
#..............................................................
# bring in mtdt for names
#..............................................................
clsts <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("hv001")) %>%
  dplyr::filter(!duplicated(.))
#..............................................................
# bring in param master list
#.............................................................
mastermap <- readRDS("analyses/cluster_inbreeding/smpls_coione/paramset/mastermap.RDS") %>%
  dplyr::select(c("1", "m", "m_learningrate", "f_learningrate", "inputpath", "parampath")) %>% # note, fstart are same across other params
  dplyr::rename(param_set = parampath,
                fstart = 1,
                mstart = m) %>%
  dplyr::mutate(
    spacetype = sub("_gens.RDS", "", sub("data/derived_data/coione_clst_inbreeding_dat/", "", inputpath))
  ) %>%
  dplyr::select(-c("inputpath"))

#..............................................................
# read in data from slurm scr
#..............................................................
filepaths <- list.files("results/clust_inbd_results/smpls_coione/Find_grad_descent_results/",
                        pattern = ".RDS", full.names = T)

read_cost_results <- function(path, clstnames){
  dat <- readRDS(path)
  # just pull out costs and see if cost was at end (or if learning rate was potentially too large) -- monotonic dec
  mincost <- min(dat$cost)
  monoton_cost_dec <- !any(diff(dat$cost) > 0)

  # make out
  out <- tibble::tibble(
    param_set = sub(".RDS", "", basename(path)),
    mincost = mincost,
    monoton_cost_dec = monoton_cost_dec
  )
  return(out)
}


mastermap_rets <- purrr::map(filepaths, read_cost_results) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(mastermap, ., by = "param_set")


#............................................................
# Process Results to extract Min Cost
#...........................................................
mastermap_mincost <- mastermap_rets %>%
  dplyr::group_by(spacetype) %>%
  dplyr::filter( mincost == min(mincost) )

#............................................................
# pull out min results
#...........................................................
clust_inb <- mastermap_mincost %>%
  dplyr::mutate(param_set = paste0("results/clust_inbd_results/smpls_coione/Find_grad_descent_results/", param_set, ".RDS"),
                param_set = purrr::map(param_set, readRDS),
                cost = purrr::map(param_set, "cost"))
#......................
# check by eye that Grad Descent looks reasonable
#......................
plot(clust_inb$cost[[1]])
plot(clust_inb$cost[[2]])

#......................
# process final Inbreeding coeffs
#......................
clust_inb$inbreed_ests <- purrr::map(clust_inb$param_set, function(prmst) {
                                                          tibble::tibble(param = c(prmst$deme_key$Deme, "m"),
                                                                         est = c(prmst$Final_Fis, prmst$Final_m))
                                                           })


#..............................................................
# out
#..............................................................
dir.create("results/clust_inbd_results/smpls_coione/min_cost_inbreedingresults/", recursive = TRUE)
write_rds(x = clust_inb,
          path = "results/clust_inbd_results/smpls_coione/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS")








