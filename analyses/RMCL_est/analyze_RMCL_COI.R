## .................................................................................
## Purpose: Look at COI estimation from RMCL from big barcode manuscript
##
##
## Notes:
## .................................................................................
library(tidyverse)
#............................................................
##### Data Wrangle ####
#...........................................................
coi <- readRDS("data/raw_data/RMCL_results/non_summariased_cois.rds") %>%
  dplyr::filter(region_denom == 1) %>% # just subset to DRc
  dplyr::filter(region == "DRC") %>%
  dplyr::filter(gt == 0.1) %>% # same as in Big Barcode Manuscript
  dplyr::mutate(barcode = ifelse(nchar(name) == 9, paste(strsplit(name, split = "")[[1]][5:nchar(name)], collapse = ""),
                                 ifelse(nchar(name) == 8, paste(strsplit(name, split = "")[[1]][4:nchar(name)], collapse = ""),
                                        ifelse(nchar(name) == 6, paste(strsplit(name, split = "")[[1]][2:nchar(name)], collapse = ""),
                                               name))))

# metadata
drcsmpls <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("barcode", "hv001", "longnum", "latnum"))

# bring together
coi_ge <- dplyr::left_join(coi, drcsmpls, by = "barcode")


#............................................................
##### Raw Maps ####
#...........................................................
# drc provs
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))
drccites <- readr::read_csv("data/map_bases/DRC_city_coordinates.csv") %>%
  dplyr::filter(population > 350000)

# by sample
ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  geom_jitter(data = coi_ge,
              aes(x = longnum, y = latnum, color = mean),
              width = 0.5, height = 0.5, alpha = 0.8) +
  scale_color_viridis_c("COI by Smpl") +
  xlab("") + ylab("") +
  theme_minimal() +
  coord_sf(datum = NA)

# by cluster mean
rmcl_summ <- coi_ge %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(longnum = mean(longnum),
                   latnum = mean(latnum),
                   stdev = sd(mean),
                   n = dplyr::n(),
                   mean = mean(mean))
rmcl_summ_plotObj <- rmcl_summ %>%
  ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  geom_point(aes(x = longnum, y = latnum, color = mean), alpha = 0.8) +
  scale_color_viridis_c("Mean COI by Clst") +
  geom_point(data = drccites,
             aes(x = longnum, y=latnum), color = "#969696") +
  ggrepel::geom_text_repel(data = drccites, aes(label = city, x = longnum, y=latnum),
                           hjust = 0.5, vjust = 0.5, nudge_y = 0.3,
                           fontface = "bold", color = "#969696") +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_sf(datum = NA)


# look at incidence
pfincidence <- raster::raster("data/raw_data/MAPrasters/getRaster/2019_Global_Pf_Incidence.201310_.18_40_8_2020_09_18.tiff")
pfincidence <- raster::mask(pfincidence, DRCprov)
pfincidence_plotObj <- ggplot() +
  ggspatial::layer_spatial(data = pfincidence, aes(fill = stat(band1))) +
  geom_sf(data = DRCprov, color = "#737373", fill = NA, size = 0.05) +
  scale_fill_distiller("Pf Incidence", type = "div", palette = "RdYlBu", na.value = NA) +
  geom_point(data = drccites,
             aes(x = longnum, y=latnum), color = "#969696") +
  ggrepel::geom_text_repel(data = drccites, aes(label = city, x = longnum, y=latnum),
                           hjust = 0.5, vjust = 0.5, nudge_y = 0.3,
                           fontface = "bold", color = "#969696") +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_sf(datum = NA)

# come together
cowplot::plot_grid(rmcl_summ_plotObj, pfincidence_plotObj, nrow = 1,
                   rel_widths = c(0.7, 1.1))


#............................................................
# IDW means
#...........................................................
library(gstat)
library(raster)
# get internal pieces needed for interpolation
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14))
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("longnum","latnum")

# need this for bounding box
DRC <- sf::as_Spatial(osmdata::getbb("Democratic Republic of the Congo",
                                     featuretype = "country",
                                     format_out = 'sf_polygon'))
drcrstr <- raster::rasterFromXYZ(cbind(grid.pred[,1],
                                       grid.pred[,2],
                                       NA),
                                 crs="+proj=longlat +datum=WGS84 +no_defs")


#......................
# quick interpolate
#......................
rmcl_summ.sf <- sf::st_as_sf(x = rmcl_summ,
                             coords = c("longnum", "latnum"),
                             crs = "+proj=longlat +datum=WGS84 +no_def")

rmcl_summ.sp <- sf::as_Spatial(rmcl_summ.sf)

gs <- gstat::gstat(formula = mean ~ 1,
                   locations = rmcl_summ.sp)
v <- variogram(gs)
plot(v)

## [inverse distance weighted interpolation]
ret <- raster::interpolate(drcrstr, idwmod, xyNames = c("longnum", "latnum"))
ret <- raster::mask(ret, DRC)
plot(ret)

#............................................................
# compare rmcl to inbreeding
#...........................................................
ge <- readRDS("data/derived_data/spacemips_GE.rds") %>%
  dplyr::mutate(hv001 = as.character(hv001))
clst_inbd <- readRDS("results/clust_inbd_results/min_cost_inbreedingresults/min_cost_inbreedingresults.RDS") %>%
  dplyr::filter(spacetype == "gcdist") %>%
  dplyr::select(c("spacetype", "inbreed_ests")) %>%
  tidyr::unnest(cols = inbreed_ests) %>%
  dplyr::mutate(hv001 = param) %>%
  dplyr::left_join(., ge, by = "hv001")

# clst inbd results
clst_inbd_plotObj <- clst_inbd %>%
  ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  geom_point(aes(x = longnum, y = latnum, color = est), alpha = 0.8) +
  scale_color_viridis_c("DISC") +
  geom_point(data = drccites,
             aes(x = longnum, y=latnum), color = "#969696") +
  ggrepel::geom_text_repel(data = drccites, aes(label = city, x = longnum, y=latnum),
                           hjust = 0.5, vjust = 0.5, nudge_y = 0.3,
                           fontface = "bold", color = "#969696") +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_sf(datum = NA)

# come together
cowplot::plot_grid(clst_inbd_plotObj, rmcl_summ_plotObj, nrow = 1)

#......................
# quick look
#......................
# cluster level
rmcl_summ %>%
  dplyr::mutate(hv001 = as.character(hv001)) %>%
  dplyr::full_join(., clst_inbd, by = "hv001") %>%
  ggplot() +
  geom_point(aes(x = mean, y = est, size = n), alpha = 0.5) +
  geom_smooth(aes(x = mean, y = est, weight = n), method = loess, show.legend = FALSE) +
  xlab("Mean Cluster COI") + ylab("DISC")

# sample level
coi_ge %>%
  dplyr::mutate(hv001 = as.character(hv001)) %>%
  dplyr::left_join(., clst_inbd, by = "hv001") %>%
  ggplot() +
  geom_point(aes(x = mean, y = est)) +
  xlab("Smpl COI") + ylab("DISC")

