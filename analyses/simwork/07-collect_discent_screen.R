## .................................................................................
## Purpose: Collect simulation results
##
## Notes:
## .................................................................................
library(tidyverse)
source("R/themes.R")

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
# pull out min costs
#...........................................................
out <- sim_disc_map %>%
  dplyr::group_by(q) %>%
  dplyr::filter(mincost == min(mincost, na.rm = T))

# save
dir.create("results/simulated_min_cost_inbreedingresults/", recursive = T)
saveRDS(out, "results/simulated_min_cost_inbreedingresults/sim_min_cost_inbreedingresults.RDS")


#............................................................
# check for convergence
#...........................................................
p1 <- out %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::mutate(costdiff = c(diff(cost), NA)) %>%
  dplyr::group_by(q) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                q = factor(q, levels = c("mq_bad", "mq_good", "ne", "coi"),
                                 labels = c("MQ-B", "MQ-G", "Ne", "COI"))) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = costdiff)) +
  facet_wrap(~q, scales = "free") +
  xlab("Iteration") +
  ylab("Cost Difference") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

p2 <- out %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::mutate(costdiff = c(diff(cost), NA)) %>%
  dplyr::group_by(q) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                q = factor(q, levels = c("mq_bad", "mq_good", "ne", "coi"),
                                 labels = c("MQ-B", "MQ-G", "Ne", "COI"))) %>%
  dplyr::filter(iteration != 1) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = costdiff)) +
  facet_wrap(~q, scales = "free") +
  xlab("Iteration") +
  ylab("Cost Difference") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

p3 <- out %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::group_by(q) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                q = factor(q, levels = c("mq_bad", "mq_good", "ne", "coi"),
                                 labels = c("MQ-B", "MQ-G", "Ne", "COI"))) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = cost)) +
  facet_wrap(~q, scales = "free") +
  xlab("Iteration") +
  ylab("Cost") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

p4 <- out %>%
  tidyr::unnest(cols = "cost") %>%
  dplyr::group_by(q) %>%
  dplyr::mutate(iteration = 1:dplyr::n(),
                q = factor(q, levels = c("mq_bad", "mq_good", "ne", "coi"),
                           labels = c("MQ-B", "MQ-G", "Ne", "COI"))) %>%
  dplyr::filter(iteration != 1) %>%
  ggplot() +
  geom_point(aes(x = iteration, y = cost)) +
  facet_wrap(~q, scales = "free") +
  xlab("Iteration") +
  ylab("Cost") +
  plot_theme +
  theme(panel.grid = element_line(color = "#bdbdbd", size = 0.1),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

p1
p2
p3
p4
