## .................................................................................
## Purpose:
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
    coord_cartesian(xlim = c(0.5,1), ylim = c(0,0.15)) +
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
# calculate and bring in distance
#...........................................................
locatdist <- dist(locats[,1:2])
locatdist <- as.matrix(locatdist)
colnames(locatdist) <- locats$deme
rownames(locatdist) <- locats$deme
locatdist_long <- locatdist %>%
  cbind.data.frame(deme = rownames(locatdist), locatdist) %>%
  tidyr::pivot_longer(., cols = -c("deme"), names_to = "deme2", values_to = "geodist") %>%
  dplyr::rename(deme1 = deme)
locatdist_long <- expand_distance_matrix(locatdist_long) %>%
  dplyr::mutate(deme1 = as.numeric(deme1),
                deme2 = as.numeric(deme2))

#......................
# join in gen and geodist
#......................
fix_gen_geo <- function(gendat, geodat) {
  dplyr::left_join(gendat, geodat, by = c("deme1", "deme2")) %>%
    dplyr::rename(gendist = simIBD,
                  locat1 = deme1,
                  locat2 = deme2) %>%
    dplyr::select(c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))
}

retmap$gengeodat <- purrr::map(retmap$data, fix_gen_geo, geodat = locatdist_long)








#............................................................
#
#...........................................................
start_params <- rep(0.1, 350)
names(start_params) <- 1:350
start_params <- c(start_params, "m" = 1e-14)
ret <- discent::deme_inbreeding_spcoef(K_gendist_geodist = retmap$gengeodat[[3]],
                                       start_params = start_params,
                                       m_lowerbound = 1e-25,
                                       m_upperbound = 5e-4,
                                       f_learningrate = 1e-8,
                                       m_learningrate = 1e-12,
                                       full_matrix = F,
                                       steps = 1e4,
                                       report_progress = TRUE)


#......................
# f values
#......................
fvals <- tibble::tibble(deme = ret$deme_key$Deme,
                        Finbd = ret$Final_Fis) %>%
  dplyr::left_join(., locats)


#............................................................
# various checks
#...........................................................
pointplot <- ggplot() +
  geom_point(data = fvals, aes(x = longnum, y = latnum, color = Finbd), size = 3, alpha = 0.5) +
  scale_color_viridis_c()


fvals <- fvals %>%
  dplyr::mutate(Findb_logit = logit(Finbd))

# prevmap fit
prevmapfit <- PrevMap::linear.model.MLE(formula = as.formula("Findb_logit ~ 1"),
                                 coords = as.formula("~longnum + latnum"),
                                 data = fvals,
                                 start.cov.pars = c(1, 1),
                                 kappa = 0.5)

#............................................................
# prevmap raster plot
#...........................................................
grid.pred <- expand.grid(c(1:300), c(1:300))
ret <- PrevMap::spatial.pred.linear.MLE(prevmapfit,
                                        grid.pred = grid.pred,
                                        scale.predictions = "prevalence",
                                        n.sim.prev = 1e2,
                                        standard.errors = T)
# make raster
ret.rstr <- raster::rasterFromXYZ(cbind(ret$grid.pred[,1],
                                        ret$grid.pred[,2],
                                        ret$prevalence$predictions),
                                  crs="+proj=longlat +datum=WGS84")

# plot out
ret.smrstr.m.plot <- ggplot() +
  ggspatial::layer_spatial(data = ret.rstr,
                           aes(fill = stat(band1)),
                           alpha = 0.8) +
  scale_fill_viridis_b("DISC", na.value = NA) +
  theme_void() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0),"cm"))


ret.smrstr.m.plot
pointplot
