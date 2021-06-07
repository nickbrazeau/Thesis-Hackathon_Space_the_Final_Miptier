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
  dplyr::filter(mincost == min(mincost, na.rm = T))


dir.create("results/cluster_inbreed_ests/min_cost_inbreedingresults/", recursive = T)
saveRDS(out, "results/cluster_inbreed_ests/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS")


#............................................................
# check for convergence
#...........................................................
p1 <- all_disc %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::mutate(costdiff = c(diff(cost), NA)) %>%
  dplyr::group_by(distlvl) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                distlvl = factor(distlvl, levels = c("gcdist", "roaddist", "airdist"),
                                 labels = c("GC Dist.", "Road Dist.", "Airport Dist."))) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = costdiff)) +
  facet_wrap(~distlvl, scales = "free") +
  xlab("Iteration") +
  ylab("Cost Difference") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

p2 <- all_disc %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::mutate(costdiff = c(diff(cost), NA)) %>%
  dplyr::group_by(distlvl) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                distlvl = factor(distlvl, levels = c("gcdist", "roaddist", "airdist"),
                                 labels = c("GC Dist.", "Road Dist.", "Airport Dist."))) %>%
  dplyr::filter(iteration != 1) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = costdiff)) +
  facet_wrap(~distlvl, scales = "free") +
  xlab("Iteration") +
  ylab("Cost Difference") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

p3 <- all_disc %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::mutate(costdiff = c(diff(cost), NA)) %>%
  dplyr::group_by(distlvl) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                distlvl = factor(distlvl, levels = c("gcdist", "roaddist", "airdist"),
                                 labels = c("GC Dist.", "Road Dist.", "Airport Dist."))) %>%
  dplyr::filter(iteration %in% 2000:4000) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = costdiff)) +
  facet_wrap(~distlvl, scales = "free") +
  xlab("Iteration") +
  ylab("Cost Difference") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

p4 <- all_disc %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::mutate(costdiff = c(diff(cost), NA)) %>%
  dplyr::group_by(distlvl) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                distlvl = factor(distlvl, levels = c("gcdist", "roaddist", "airdist"),
                                 labels = c("GC Dist.", "Road Dist.", "Airport Dist."))) %>%
  dplyr::filter(iteration %in% 9000:1e4) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = costdiff)) +
  facet_wrap(~distlvl, scales = "free") +
  xlab("Iteration") +
  ylab("Cost Difference") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))


cowplot::plot_grid(p1, p2, p3, p4, nrow = 4, align = "v")
