## .................................................................................
## Purpose:
##
##
## Notes:
## .................................................................................
library(tidyverse)

#............................................................
# functions
#...........................................................
make_ibd_histogram <- function(ibddf) {
  mainplot <- ibddf %>%
    ggplot() +
    geom_histogram(aes(x=simIBD, y = (..count../sum(..count..))*100),
                   color = "#000000", fill = "#d9d9d9") +
    xlab("IBD") + ylab("frequency (%)") +
    theme_classic()

  insetplot <- ibddf %>%
    ggplot() +
    geom_histogram(aes(x=simIBD, y = (..count../sum(..count..))*100),
                   color = "#000000", fill = "#d9d9d9") +
    xlab("IBD") + ylab("frequency (%)") +
    theme_classic() +
    coord_cartesian(xlim = c(0.5,1)) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"))
  # out
  cowplot::ggdraw() +
    cowplot::draw_plot(mainplot, x = 0, y = 0, width = 1, height = 1, scale = 1) +
    cowplot::draw_plot(insetplot, x = 0.5, y= 0.3, width = 0.4, height = 0.4)
}


#............................................................
# read in
#...........................................................
simdat <- readRDS("data/sim_data/swf_simulations.rds")
locats <- readRDS("data/sim_data/sim_locations.rds")

#............................................................
# get ibd sim results
#...........................................................
retfiles <- list.files("results/sim_clust_results/swf_sim_ret/",
                       full.names = T)
retmap <- tibble::tibble(lvl = stringr::str_split_fixed(basename(retfiles), "_", 2)[,1],
                         simIBD = retfiles) %>%
  dplyr::mutate(simIBD = purrr::map(simIBD, readRDS)) %>%
  tidyr::unnest(cols = simIBD) %>%
  dplyr::group_by(lvl) %>%
  tidyr::nest()

#......................
# get ibd plots
#......................
retmap <- retmap %>%
  dplyr::mutate(ibd_plotObj = purrr::map(data, make_ibd_histogram))

retmap$ibd_plotObj[[1]]
retmap$ibd_plotObj[[2]]
retmap$ibd_plotObj[[3]]




sum(retmap$data[[1]]$simIBD != 0) / nrow(retmap$data[[1]])
ibd <- retmap$data[[1]][retmap$data[[1]]$simIBD != 0, ]

strngth <- ibd %>%
  dplyr::group_by(deme1) %>%
  dplyr::summarise(ibdsum = sum(simIBD)) %>%
  dplyr::rename(deme = deme1)



locats_sm <- locats %>%
  dplyr::filter(deme %in% unique(ibd$deme1)) %>%
  dplyr::left_join(., strngth)
ggplot() +
  geom_tile(data = simdat$gridmig[[1]], aes(x = longnum, y = latnum, fill = migration)) +
  geom_point(data = locats_sm, aes(x = longnum, y = latnum, color = factor(deme),
                                   size = ibdsum), alpha = 0.8, show.legend = F) +
  scale_fill_viridis_c("Migration Topology") +
  coord_fixed() +
  theme_void()
