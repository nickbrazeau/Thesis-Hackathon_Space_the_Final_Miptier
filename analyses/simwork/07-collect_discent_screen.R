## .................................................................................
## Purpose: Collect simulation results
##
## Notes:
## .................................................................................
library(tidyverse)

#............................................................
# read in simulation discent results
#...........................................................

# get paths and read in results
sim_disc_map <- tibble::tibble(paths = list.files("results/sim_cluster_inbreed_ests/",
                                            full.names = T)) %>%
  dplyr::mutate(ret = purrr::map(paths, readRDS)) %>%
  tidyr::unnest(cols = "ret") %>%
  dplyr::mutate(cost = purrr::map(discentret, "cost"),
                mincost = purrr::map_dbl(cost, function(x){x[length(x)]})) %>%
  dplyr::select(-c("paths"))

#......................
# look at cost distributions
#......................
sim_disc_map %>%
  dplyr::select(c("q", "mincost")) %>%
  ggplot() +
  geom_boxplot(aes(x = mincost)) +
  coord_flip() +
  facet_wrap(~q, scales = "free")



#............................................................
# save out min costs
#...........................................................
out <- sim_disc_map %>%
  dplyr::group_by(q) %>%
  dplyr::filter(mincost == min(mincost, na.rm = T))


dir.create("results/simulated_min_cost_inbreedingresults/", recursive = T)
saveRDS(out, "results/simulated_min_cost_inbreedingresults/sim_min_cost_inbreedingresults.RDS")


