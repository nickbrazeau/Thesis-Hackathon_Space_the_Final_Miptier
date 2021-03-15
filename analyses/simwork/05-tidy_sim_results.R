## .................................................................................
## Purpose: Pull in cluster results and tidy them up
##
##
## Notes:
## .................................................................................
library(tidyverse)
source("R/pairwise_helpers.R")

#............................................................
# functions
#...........................................................
make_ibd_histogram <- function(ibddf, filt) {
  mainplot <- ibddf %>%
    dplyr::filter(simIBD > filt) %>%
    ggplot() +
    geom_histogram(aes(x=simIBD, y = (..count../sum(..count..))*100),
                   color = "#000000", fill = "#d9d9d9") +
    xlab("IBD") + ylab("frequency (%)") +
    theme_classic()

  insetplot <- ibddf %>%
    dplyr::filter(simIBD > filt) %>%
    ggplot() +
    geom_histogram(aes(x=simIBD, y = (..count../sum(..count..))*100),
                   color = "#000000", fill = "#d9d9d9") +
    xlab("IBD") + ylab("frequency (%)") +
    theme_classic() +
    coord_cartesian(xlim = c(0.5,1), ylim = c(0,1)) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"))
  # out
  cowplot::ggdraw() +
    cowplot::draw_plot(mainplot, x = 0, y = 0, width = 1, height = 1, scale = 1) +
    cowplot::draw_plot(insetplot, x = 0.5, y= 0.3, width = 0.4, height = 0.4)
}


#............................................................
# get ibd sim results
#...........................................................
retfiles <- list.files("results/swf_sim_ret/",
                       full.names = T)
retmap <- tibble::tibble(lvl = stringr::str_split_fixed(basename(retfiles), "_", 2)[,1],
                         simIBD = retfiles) %>%
  dplyr::mutate(simIBD = purrr::map(simIBD, readRDS)) %>%
  tidyr::unnest(cols = simIBD) %>%
  dplyr::group_by(lvl) %>%
  tidyr::nest() %>%
  dplyr::ungroup()

#......................
# get ibd plots
#......................
retmap <- retmap %>%
  dplyr::mutate(ibd_plotObj_all = purrr::map(data, make_ibd_histogram, filt = -1),
                ibd_plotObj_nonzero =  purrr::map(data, make_ibd_histogram, filt = 0))

retmap$ibd_plotObj_all[[1]]
retmap$ibd_plotObj_all[[2]]
retmap$ibd_plotObj_all[[3]]

retmap$ibd_plotObj_nonzero[[1]]
retmap$ibd_plotObj_nonzero[[2]]
retmap$ibd_plotObj_nonzero[[3]]

#............................................................
# bring in distance
#   remember considering m-dst twice
#...........................................................
#......................
# regular sim datas
#......................
distmap_reg <- tibble::tibble(lvl = c("coi", "ne", "mq", "mq"),
                              q = c("coi", "ne", "mq_good", "mq_bad"),
                              path = c("data/sim_data/euclidean_geodist.rds",
                                       "data/sim_data/euclidean_geodist.rds",
                                       "data/sim_data/gridmig_geodist.rds",
                                       "data/sim_data/euclidean_geodist.rds"
                              ))  %>%
  dplyr::mutate(geodist = purrr::map(path, readRDS)) %>%
  dplyr::select(-c("path"))

#......................
# non-linear sim datas
#   need to slight liftover here
#......................
distmap_nonlinear <- tibble::tibble(lvl = c("mtn", "rift", "oppcorner"),
                                    q = c("mtn", "rift", "oppcorner"),
                                    path = c("data/sim_data/mtn_nonlinear_migration_geodist.rds",
                                             "data/sim_data/rift_nonlinear_migration_geodist.rds",
                                             "data/sim_data/oppcorner_nonlinear_migration_geodist.rds"
                                    ))  %>%
  dplyr::mutate(geodist = purrr::map(path, readRDS)) %>%
  dplyr::select(-c("path")) %>%
  dplyr::mutate(geodist =  purrr::map(geodist, function(x){
    # liftover function from a distance matrix wide to long format and rename demes
    ret <- broom::tidy(as.dist(x))
    colnames(ret) <- c("deme1", "deme2", "distval")
    # make sure we are dealing w/ numerics and not factors here for deme names
    ret <- ret %>%
      dplyr::mutate(deme1 = as.integer(as.character(deme1)),
                    deme2 = as.integer(as.character(deme2)))
    return(ret)
  }))


#......................
# bring together
#......................
retmap <- dplyr::bind_rows(distmap_reg, distmap_nonlinear) %>%
  dplyr::left_join(retmap, ., by = "lvl") %>%
  dplyr::rename(gendat = data)


#......................
# join in gen and geodist
#......................
fix_gen_geo <- function(gendat, geodist) {
  dplyr::left_join(gendat, geodist, by = c("deme1", "deme2")) %>%
    dplyr::rename(gendist = simIBD,
                  geodist = distval,
                  locat1 = deme1,
                  locat2 = deme2) %>%
    dplyr::mutate(geodist = ifelse(locat1 == locat2, 0, geodist)) %>%  # NB, when save these out in tidy, no self comparison -- need that here (was on diagonal that we used for matrix)
    dplyr::select(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))
}

retmap$gengeodat <- purrr::pmap(retmap[,c("gendat", "geodist")], fix_gen_geo)


#............................................................
# Plot realized IBD pairings
#...........................................................
#......................
# plotting function
#......................
plot_ibd_pairs <- function(gengeodat, locats) {
  # quick manips
  locats1 <- locats %>%
    dplyr::rename(locat1 = deme)
  locats2 <- locats %>%
    dplyr::rename(locat2 = deme)
  # bring together
  gengeodat %>%
    dplyr::left_join(., y = locats1, by = "locat1") %>%
    dplyr::left_join(., y = locats2, by = "locat2") %>%
    dplyr::filter(gendist > 0 ) %>%
    ggplot() +
    geom_segment(aes(x = longnum.x,
                     y = latnum.x,
                     xend = longnum.y,
                     yend = latnum.y,
                     color = gendist),
                 alpha = 0.5) +
    scale_color_viridis_c("Sim. IBD") +
    theme_void()

}

# read in locats
locats <- readRDS("data/sim_data/lattice_model.rds") %>%
  dplyr::select(-c("migration"))

# run function
retmap <- retmap %>%
  dplyr::mutate(pairPlotObj = purrr::map(gengeodat, plot_ibd_pairs, locats = locats))

retmap$pairPlotObj[[1]]
retmap$pairPlotObj[[2]]
retmap$pairPlotObj[[3]]
retmap$pairPlotObj[[4]]
retmap$pairPlotObj[[5]]
retmap$pairPlotObj[[6]]
retmap$pairPlotObj[[7]]


#......................
# save out
#......................
saveRDS(retmap, "data/sim_data/sim_gengeodat.rds")

