## .................................................................................
## Purpose: Get MinCost of DRC DISC fits
##
##
## Notes:
## .................................................................................
library(tidyverse)

#............................................................
# read in results
#...........................................................
ret_read_wrapper <- function(path) {
  readRDS(path) %>%
    dplyr::mutate(datalvl = stringr::str_extract(inputpath, "allsmpls|coione"),
                  distlvl = stringr::str_extract(inputpath, "gcdist|roaddist|airdist"))
}

# get paths and read in results
retmap <- tibble::tibble(paths = list.files("results/cluster_inbreed_ests/",
                                            full.names = T)) %>%
  dplyr::mutate(ret = purrr::map(paths, ret_read_wrapper)) %>%
  tidyr::unnest(cols = "ret") %>%
  dplyr::mutate(cost = purrr::map(discentret, "cost"),
                mincost = purrr::map_dbl(cost, function(x){x[length(x)]})) %>%
  dplyr::select(-c("paths"))

#......................
# look at cost distributions
#......................
retmap %>%
  dplyr::select(c("datalvl", "distlvl", "mincost")) %>%
  ggplot() +
  geom_boxplot(aes(x = mincost)) +
  coord_flip() +
  facet_wrap(datalvl ~ distlvl, scales = "free")



#............................................................
# save out min costs
#...........................................................
out <- retmap %>%
  dplyr::group_by(datalvl, distlvl) %>%
  dplyr::filter(mincost == min(mincost))

# drop extra
out <- out %>%
  dplyr::select(c("datalvl", "distlvl", "discentret", "mincost"))


dir.create("results/cluster_inbreed_ests/min_cost_inbreedingresults/", recursive = T)
saveRDS(out, "results/cluster_inbreed_ests/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS")
