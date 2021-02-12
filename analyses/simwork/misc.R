source("R/basics.R")



#............................................................
#
#...........................................................
start <- Sys.time()
start_params <- rep(0.1, 350)
names(start_params) <- 1:350
start_params <- c(start_params, "m" = 0.1)
retdisc <- discent::deme_inbreeding_spcoef(K_gendist_geodist = retmap$gengeodat[[1]],
                                           start_params = start_params,
                                           m_lowerbound = 1e-10,
                                           m_upperbound = 1,
                                           f_learningrate = 1e-5,
                                           m_learningrate = 1e-10,
                                           full_matrix = F,
                                           steps = 5e4,
                                           report_progress = TRUE)

Sys.time() - start
#......................
# f values
#......................
fvals <- tibble::tibble(deme = retdisc$deme_key$Deme,
                        Finbd = retdisc$Final_Fis) %>%
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
