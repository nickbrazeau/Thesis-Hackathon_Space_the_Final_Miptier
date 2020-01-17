source("R/pairwise_helpers.R")
source("R/themes.R")
mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("name", "barcode", "hv001", "longnum", "latnum"))

ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::select(c("hv001", "urban_rura")) %>%
  dplyr::mutate(urban_rura_ext = ifelse(urban_rura == "R", 10000, 2000))

mtdt <- dplyr::left_join(mtdt, ge, by = "hv001")
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD.long.mtdt.rds")

#..............................................................
# urbanicity
#..............................................................
urban <- raster::raster("data/derived_data/urbanicity_raster/urbanicity.grd")
urbanmean <- rep(NA, nrow(mtdt))
for (i in 1:nrow(mtdt)){
  urbanmean[i] <- raster::extract(x = urban,
                                  y = sf::as_Spatial(mtdt$geometry[i]),
                                  buffer = mtdt$urban_rura_ext[i],
                                  fun = mean,
                                  na.rm = T)
}
mtdt$urban <- urbanmean

#..............................................................
# all between
#..............................................................

ibD.long.btwn <- ibD %>%
  dplyr::filter(hv001.x != hv001.y) %>%
  expand_distance_matrix(.) %>%
  dplyr::select(c("smpl1", "smpl2", "malecotf", "hv001.x", "hv001.y"))
node.edge.dens.weighted.btwn <- ibD.long.btwn %>%
  dplyr::group_by(smpl1) %>%
  dplyr::summarise(
    n = n(),
    F_wi = sum(malecotf)) %>% # \sum edge * malecotf (NB each edge is 1)
  dplyr::rename(name = smpl1)


node.edge.dens.weighted.btwn <- dplyr::left_join(node.edge.dens.weighted.btwn,
                                                 mtdt, by = "name")

plotObj.btwn <- node.edge.dens.weighted.btwn %>%
  ggplot() +
  geom_point(aes(x = urban, y = F_wi)) +
  plot_theme +
  xlab("Urbanicity Factor Score") +
  ylab("Weighted Edge Density")

node.edge.dens.weighted.btwn.dcortest <-
  energy::dcor.test(node.edge.dens.weighted.btwn$F_wi,
                    node.edge.dens.weighted.btwn$urban,
                    R = 1e4)




#..............................................................
# meiotics between
#..............................................................
ibD.long.meioticbtwn <- ibD %>%
  dplyr::filter(malecotf > 0.5) %>%
  dplyr::filter(hv001.x != hv001.y) %>%
  expand_distance_matrix(.) %>%
  dplyr::select(c("smpl1", "smpl2", "malecotf", "hv001.x", "hv001.y"))
node.edge.dens.weighted.meioticbtwn <- ibD.long.meioticbtwn %>%
  dplyr::group_by(smpl1) %>%
  dplyr::summarise(
    n = n(),
    F_wi = sum(malecotf)) %>% # \sum edge * malecotf (NB each edge is 1)
  dplyr::rename(name = smpl1)


node.edge.dens.weighted.meioticbtwn <- dplyr::left_join(node.edge.dens.weighted.meioticbtwn,
                                                         mtdt, by = "name")

plotObj.meioticbtwn <- node.edge.dens.weighted.meioticbtwn %>%
  ggplot() +
  geom_point(aes(x = urban, y = F_wi)) +
  plot_theme +
  xlab("Urbanicity Factor Score") +
  ylab("Weighted Edge Density")

node.edge.dens.weighted.meioticbtwn.dcortest <-
  energy::dcor.test(node.edge.dens.weighted.meioticbtwn$F_wi,
                    node.edge.dens.weighted.meioticbtwn$urban,
                    R = 1e4)


#..............................................................
# all
#..............................................................
ibD.long.all <- ibD %>%
  expand_distance_matrix(.) %>%
  dplyr::select(c("smpl1", "smpl2", "malecotf", "hv001.x", "hv001.y"))
node.edge.dens.weighted.all <- ibD.long.all %>%
  dplyr::group_by(smpl1) %>%
  dplyr::summarise(
    n = n(),
    F_wi = sum(malecotf)) %>% # \sum edge * malecotf (NB each edge is 1)
  dplyr::rename(name = smpl1)


node.edge.dens.weighted.all <- dplyr::left_join(node.edge.dens.weighted.all,
                                                mtdt, by = "name")

plotObj.all <- node.edge.dens.weighted.all %>%
  ggplot() +
  geom_point(aes(x = urban, y = F_wi)) +
  plot_theme +
  xlab("Urbanicity Factor Score") +
  ylab("Weighted Edge Density")

node.edge.dens.weighted.all.dcortest <-
  energy::dcor.test(node.edge.dens.weighted.all$F_wi,
                    node.edge.dens.weighted.all$urban,
                    R = 1e4)



#..............................................................
# meiotics all
#..............................................................
ibD.long.meioticall <- ibD %>%
  dplyr::filter(malecotf > 0.5) %>%
  expand_distance_matrix(.) %>%
  dplyr::select(c("smpl1", "smpl2", "malecotf", "hv001.x", "hv001.y"))
node.edge.dens.weighted.meioticall <- ibD.long.meioticall %>%
  dplyr::group_by(smpl1) %>%
  dplyr::summarise(
    n = n(),
    F_wi = sum(malecotf)) %>% # \sum edge * malecotf (NB each edge is 1)
  dplyr::rename(name = smpl1)


node.edge.dens.weighted.meioticall <- dplyr::left_join(node.edge.dens.weighted.meioticall,
                                                         mtdt, by = "name")

plotObj.meioticall <- node.edge.dens.weighted.meioticall %>%
  ggplot() +
  geom_point(aes(x = urban, y = F_wi)) +
  plot_theme +
  xlab("Urbanicity Factor Score") +
  ylab("Weighted Edge Density")

node.edge.dens.weighted.meioticall.dcortest <-
  energy::dcor.test(node.edge.dens.weighted.meioticall$F_wi,
                    node.edge.dens.weighted.meioticall$urban,
                    R = 1e4)




#..............................................................
# Figure out
#..............................................................
jpeg("~/Documents/GitHub/Space_the_Final_Miptier/results/figures/edge_dens_urban.jpg",
     height = 8, width = 8, units = "in", res = 500)
cowplot::plot_grid(plotObj.all,
                   plotObj.meioticall,
                   plotObj.btwn,
                   plotObj.meioticbtwn,
                   align = "h", nrow = 2, labels = c("(A)", "(B)", "(C)", "(D)"))
graphics.off()

jpeg("~/Documents/GitHub/Space_the_Final_Miptier/results/figures/edge_dens_urban.jpg",
     height = 8, width = 11, units = "in", res = 500)
cowplot::plot_grid(plotObj.meioticall,
                   plotObj.meioticbtwn,
                   align = "h", nrow = 1, labels = c("(A)", "(B)"))
graphics.off()

#..............................................................
# Table out
#..............................................................
dattests <- ls()[grepl("dcortest", ls())]
ret <- tibble::tibble(
  lvl = dattests
  ) %>%
  dplyr::mutate(
    lvl = stringr::str_split_fixed(lvl, "\\.", n=6)[,5],
    dat = lapply(dattests, get)
  )

ret$pvalue <- purrr::map(ret$dat, "p.value")
ret$dCor <- purrr::map(ret$dat, "statistic")

out <- ret %>%
  dplyr::select(-c("dat")) %>%
  tidyr::unnest(cols = c("pvalue", "dCor")) %>%
  dplyr::mutate_if(is.numeric, round, 4)

