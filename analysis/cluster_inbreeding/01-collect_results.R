#########################################################################
# Purpose: Collect the results from the cluster inbreeding calculation
#
# Author: Nicholas F. Brazeau
#
# Date: April 08 2020
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
mastermap.lg <- readRDS("data/derived_data/clst_inbreeding_dat/paramset/mastermap.RDS")
mastermap <- mastermap.lg %>%
  dplyr::select(c("1", "m", "learningrate", "inputpath", "parampath")) %>% # note, fstart are same across other params
  dplyr::rename(param_set = parampath,
                fstart = 1,
                mstart = m) %>%
  dplyr::mutate(
    spacetype = sub("_gens.RDS", "", sub("data/derived_data/clst_inbreeding_dat/", "", inputpath))
  ) %>%
  dplyr::select(-c("inputpath"))

#..............................................................
# read in data from slurm scr
#..............................................................
# NB cluster names are sorted for Fij results sort(clsts$hv001)
filepaths <- list.files("data/derived_data/clst_inbreeding_dat/clust_results/",
                        pattern = ".RDS", full.names = T)

read_cost_results <- function(path, clstnames){
  dat <- readRDS(path)
  # just pull out costs and see if cost was at end (or if learning rate was too large) -- monotonic dec
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


cost_rets <- purrr::map(filepaths, read_cost_results) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(mastermap, ., by = "param_set")


# save out for later use
saveRDS(clst_rets, "data/derived_data/clst_inbreeding_dat/clust_inbd_results/collective_clust_inbd_results.RDS")
#............................................................
# Process Results to see if Convergence was Reached
#...........................................................
clst_rets <- readRDS("data/derived_data/clst_inbreeding_dat/clust_inbd_results/collective_clust_inbd_results.RDS")

# going to make a seperate data frame for looking at convergence for a given param set
conv_filter_map <- clst_rets %>%
  dplyr::select(c("param_set", "conv_stats")) %>%
  dplyr::left_join(., y = mastermap, by = "param_set")

check_conv <- function(conv_stats, conv_alpha, learningrate) {
  mean(sapply(conv_stats, function(x, y) {x <= y}, y = -conv_alpha*learningrate))
}

conv_filter_map <- conv_filter_map %>%
  dplyr::mutate(conv_pass = purrr::pmap_dbl(.l = conv_filter_map[,c("conv_stats", "learningrate")],
                                            .f = check_conv,
                                            conv_alpha = 1e-5)) %>%
  dplyr::select(c("param_set", "conv_pass", "spacetype"))

# now attach to master map and fitler
mastermap.results <- dplyr::left_join(clst_rets, conv_filter_map, by = "param_set") %>%
  dplyr::filter(conv_pass >= 0.99) %>%
  dplyr::select(-c("conv_stats")) %>%
  dplyr::select(c("param_set", "spacetype", dplyr::everything()))

#..............................................................
#out
#..............................................................
write_rds(x = mastermap.results,
          path = "data/derived_data/clst_inbreeding_dat/clst_inbreeding_convfilt_results.RDS")


#..............................................................
# play
#..............................................................

mastermap.results %>%
  dplyr::select(-c("conv_pass")) %>%
  dplyr::filter(m < 0.99) %>%
  dplyr::select(-c("m")) %>%
  tidyr::gather(., key = "param", value = "est", 3:ncol(.)) %>%
  dplyr::group_by(spacetype, param) %>%
  dplyr::summarise(
    n = dplyr::n(),
    accepted = n/(nrow(mastermap.lg)/3),
    min = min(est),
    LCI = quantile(est, 0.025),
    median = median(est),
    mean = mean(est),
    UCI = quantile(est, 0.975),
    max = max(est)
  ) %>%
  ggplot() +
  geom_pointrange(aes(x = param, y = mean, ymin = LCI, ymax = UCI)) +
  facet_wrap(~ spacetype, nrow = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  )


















# sanity
