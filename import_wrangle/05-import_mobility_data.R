#########################################################################
# Purpose:
#
# Author: Nicholas F. Brazeau
#
# Data Source: World Pop Internal Flows https://www.worldpop.org/geodata/summary?id=1281
#       Methods: https://www.nature.com/articles/sdata201666
#                https://academic.oup.com/migration/article/3/1/89/2413406
#########################################################################
library(tidyverse)
library(sf)
source("R/basics.R")
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))

#............................................................
# read in data
#...........................................................
mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("name", "barcode", "hv001", "longnum", "latnum"))

flows <- readr::read_csv("data/raw_data/worldpop/COD_5yrs_InternalMigFlows_2010/COD_5yrs_InternalMigFlows_2010.csv")
centroids <- sf::read_sf("data/raw_data/worldpop/COD_5yrs_InternalMigFlows_2010/COD_AdminUnit_Centroids/COD_AdminUnit_Centroids.shp")
centroids <- sf::st_transform(centroids, sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs"))

#......................
# voronoi tesselations from migration data
#......................
vr = sf::st_sfc(sf::st_voronoi(centroids, envelope = DRCprov))


#......................
# cities for plotting/sanity check
#......................
drccities <- readr::read_csv("data/map_bases/DRC_city_coordinates.csv") %>%
  dplyr::filter(population >= 5e4)
# liftover for drc cities
drccities <- drccities %>%
  dplyr::mutate(pop_fact = cut(population, breaks = c(50e3, 100e3, 250e3, 500e3, 1e6, Inf), right = F),
                pop_fact = factor(pop_fact, labels = c("50-100", "100-250", "250-500", "500-1,000", ">1,000")))


airports <- readr::read_csv("data/raw_data/flight_data/hotosm_cd-airports.csv") %>%
  dplyr::filter(type %in% c("large_airport", "medium_airport")) %>%
  dplyr::select(c("name", "longitude_deg", "latitude_deg"))
airports <- sf::st_as_sf(airports, coords = c("longitude_deg", "latitude_deg"), crs = 4326)


drcadm <- raster::getData(name = "GADM", country = "COD", level = 1, path = tempdir())

flows.sub <- flows %>%
  dplyr::mutate(prdmigscaled = my.scale(PrdMIG)) %>%
  dplyr::filter(prdmigscaled > 0)
ggplot() +
  geom_sf(data = sf::st_as_sf(drcadm)) +
  #geom_point(data = mtdt, aes(x = longnum, y = latnum)) +
  geom_curve(data = flows.sub, alpha = 0.5, size = 0.8,
             aes(x = LONFR, y = LATFR,
                 xend = LONTO, yend = LATTO,
                 color = prdmigscaled)) +
  geom_sf(data = centroids, color = "red", size = 3, alpha = 0.5) +
  ggrepel::geom_text_repel(data = drccities, aes(label = city, x = longnum, y=latnum),
                            hjust = 0.5, vjust = 0.5, nudge_y = 0.8, fontface = "bold") +
  geom_point(data = drccities, aes(x = longnum, y=latnum), shape = 21, stroke = 0.3, fill = NA) +
  scale_color_viridis_c()





#............................................................
#
#...........................................................
library(tidyverse)
source("R/pairwise_helpers.R")

flows <- readr::read_csv("data/raw_data/worldpop/COD_5yrs_InternalMigFlows_2010/COD_5yrs_InternalMigFlows_2010.csv") %>%
  magrittr::set_colnames(tolower(colnames(.)))

nodes <- unique(c(flows$NODEI, flows$NODEJ))
lftovr <- data.frame(t(combn(nodes, 2))) %>%
  magrittr::set_colnames(c("nodei", "nodej")) %>%
  dplyr::mutate(id = 1:nrow(.))
lftovr <- expand_distance_matrix(lftovr)

# get diff
flows.diff <- flows %>%
  dplyr::left_join(., lftovr)
flows.diff <- split(flows.diff, factor(flows.diff$id))
flows.diff <- lapply(flows.diff, function(x) {
  # ij - ji
  x <- x %>%
    dplyr::arrange(nodei)
  ret <- x[1,] %>%
    dplyr::select(-c("prdmig"))
  ret$netmig <- x$prdmig[1] - x$prdmig[2]
  return(ret) }) %>%
  dplyr::bind_rows()

# scale
flows.diff$netmig <- flows.diff$netmig/sd(flows.diff$netmig)

# plotdf
flows.diff_plotdf <- flows.diff %>%
  # note, need to determine "correct end" for arrow
  dplyr::mutate(long.start = ifelse(netmig > 0, lonfr, lonto),
                lat.start = ifelse(netmig > 0, latfr, latto),
                long.end = ifelse(netmig > 0, lonto, lonfr),
                lat.end = ifelse(netmig > 0, latto, latfr),
                from = ifelse(netmig > 0, nodei, nodej),
                to = ifelse(netmig > 0, nodej, nodei))

ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  geom_curve(data = flows.diff_plotdf,
             aes(x = long.start, y = lat.start,
                 xend = long.end, yend = lat.end,
                 color = netmig),
             arrow = arrow(length = unit(0.02, "npc")),
             alpha = 0.5) +
  geom_sf(data = centroids,
             show.legend = F, shape = 21,
             fill = "#ffffff", color = "#000000" ) +
  scale_color_viridis_c("Net Migration", option="plasma", direction = 1)


# gradient
flows.diff.grad_plotdf <- flows.diff_plotdf %>%
  dplyr::filter(netmig > 0.1 | netmig < -0.1)

ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  geom_curve(data = flows.diff.grad_plotdf,
             aes(x = long.start, y = lat.start,
                 xend = long.end, yend = lat.end,
                 color = netmig),
             arrow = arrow(length = unit(0.02, "npc")),
             alpha = 0.5) +
  geom_sf(data = centroids,
          show.legend = F, shape = 21,
          fill = "#ffffff", color = "#000000" ) +
  scale_color_viridis_c("Net Migration", option="plasma", direction = 1)


centroids <- centroids %>%
  dplyr::rename(nodei = IPUMSID)

flows.diff.simple_plotdf <- flows.diff_plotdf %>%
  dplyr::select(c("nodei", "nodej", "netmig")) %>%
  expand_distance_matrix(.) %>%
  dplyr::group_by(nodei) %>%
  dplyr::summarise(meannetmig = mean(netmig),
                   sdnetmig = sd(netmig)) %>%
  dplyr::left_join(., centroids, by = "nodei")


ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  geom_point(data = flows.diff.simple_plotdf, aes(x = POINT_X, y = POINT_Y, color = meannetmig, size = sdnetmig), alpha = 0.5) +
  scale_color_viridis_c("Net Migration", option="plasma", direction = 1)
