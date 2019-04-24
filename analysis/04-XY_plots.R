#----------------------------------------------------------------------------------------------------
# Purpose of this script is to plot X-Y of space and genetic distance and then look at potential covariates
# that shape these relationships
#----------------------------------------------------------------------------------------------------
# libraries and imports
gdrive <- tcltk::tk_choose.dir()
library(tidyverse)

# grab data and example couple of plots
drcmips <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_gen_space_epi_final.rds"))

dist_grouped <- gen_spatial_dist_plot(drcmips,select = "Barcode",spatial_bins = 25)
admin_plot_mle <- gen_spatial_dist_plot(drcmips, select = "ADM1NAME", genetic_measures = "MLE_IBD",
                                    spatial_bins = NULL, facet_y = "ADM1NAME_1")

cowplot::save_plot(paste0(gdrive, "/plots/ibd.png"), admin_plot_mle + ggtitle("") + ylab("Relatedness (F)"))

# ---------------------------

# grab the data from the admin plot andbind in covariates
dat <- admin_plot_mle$data
meta <- readRDS(paste0(gdrive,"/data/derived_data/cd2013_admin_covariates.rds"))
dat <- left_join(dat, rename(meta, ADM1NAME_1 = ADM1NAME), by = "ADM1NAME_1")
dat$ADM1NAME_2 <- strsplit(as.character(dat$bin),"|",fixed=TRUE) %>% lapply(function(x) x[2]) %>% unlist

# make new covariates
dat$gen_dist <- dat$gen_dist + .Machine$double.xmin
dat$spat_dist_scale <- dat$spat_dist/max(dat$spat_dist)
dat$travel_scale <- dat$travel/max(dat$travel)
dat$altitude_scale <- dat$altitude/max(dat$altitude)

# differnece in prev
dat$prev_diff <- meta$micro_prev[match(dat$ADM1NAME_1,meta$ADM1NAME)] - meta$micro_prev[match(dat$ADM1NAME_2,meta$ADM1NAME)]

# markham seaonality index
markham_seasonality_index <- function(an_vec) {

  pos <- c(0,cumsum(lubridate::days_in_month(1:12)))
  rads <-  midpoints(pos)/365 * 2 * pi
  names(rads) <- tail(names(pos),-1)
  pos[1] <- 1

  I <- rep(0,12)
  for (i in 1:12) {
    I[i] <- sum(an_vec[pos[i]:pos[i+1]])
  }
  return(sqrt(sum(I * sin(rads))^2 + sum(I * cos(rads))^2)/sum(I))
}
vecs <- lapply(unique(dat$dhs_region),magenta:::seasonal_profile,"Democratic Republic of the Congo")
dat$msi <- unlist(lapply(vecs, markham_seasonality_index))[match(dat$dhs_region, unique(dat$dhs_region))]

# MODEL SELECTION
# ---------------------------------

# Start with what we talked about as sensible model using just the prevalence difference
pd_mod <- lme4::glmer(gen_dist ~ spat_dist_scale + prev_diff + (1+spat_dist_scale|ADM1NAME_1),
  data = dat, family=gaussian(link="log"),
  weights = length)

# build most complex sensible model
all_mod <- lme4::glmer(
  gen_dist ~ (spat_dist_scale + travel_scale + prev_diff + msi) +
    (1+spat_dist_scale|ADM1NAME_1),
  data = dat, family=gaussian(link="log"),
  weights = length)

# explore combinations
dr <- MuMIn::dredge(all_mod, m.lim = c(1,32))

# best fit and sensible m is our original pd_mod

dat <- cbind(dat, merTools::predictInterval(pd_mod,type="probability"))

# interval_plot <- admin_plot_mle +
#   geom_line(data = dat, mapping = aes(x=spat_dist,y=fit),color="red") +
#   geom_line(data = dat, mapping = aes(x=spat_dist,y=upr),color="red", linetype = "dashed") +
#   geom_line(data = dat, mapping = aes(x=spat_dist,y=lwr),color="red", linetype = "dashed")

fit_plot <- admin_plot_mle +
  geom_line(data = dat, mapping = aes(x=spat_dist,y=fit),color="red",lwd=1)

# our model suggests that a given region, r1, if r1 has higher prevalence than
# r2 then the greater the prevalence difference the less related they are. If r2 has a higher
# prevalence relationship then the inverse is found.
summary(pd_mod)

## as demonstrated with the following
newdata = expand.grid(spat_dist_scale=seq(0.01,0.99,0.01),prev_diff=seq(-0.75,0.75,0.01))
pars <- summary(pd_mod)$coefficients[,"Estimate"]
newdata$pred <- exp(pars[1] +
                      pars[2]*newdata$spat_dist_scale +
                      pars[3]*newdata$prev_diff +
                      pars[4]*newdata$spat_dist_scale*newdata$prev_diff)

breaks <- seq(min(newdata$pred), max(newdata$pred), length.out = 11)
newdata$breaks <- cut(newdata$pred, breaks,include.lowest = TRUE)
newdata$pred_mid <- as.factor(round(midpoints(breaks),3)[newdata$breaks])

cont <- ggplot(newdata, aes(x=spat_dist_scale*max(dat$spat_dist),y=prev_diff,z=pred_mid,fill=pred_mid)) +
  geom_tile(aes(fill = pred_mid)) +
  scale_fill_viridis_d(name = "Relatedness") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Microscopy Prevalence Difference") +
  xlab("Greater Circle Distance (km)")
cowplot::save_plot(paste0(gdrive, "/plots/contour.png"), cont)
