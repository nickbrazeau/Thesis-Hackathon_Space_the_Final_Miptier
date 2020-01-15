source("R/pairwise_helpers.R")
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
# adjacency matrix
#..............................................................
ibD.long.meiotic <- ibD %>%
  expand_distance_matrix(.) %>%
  dplyr::select(c("smpl1", "smpl2", "malecotf", "hv001.x", "hv001.y"))
node.edge.dens.weighted <- ibD.long.meiotic %>%
  dplyr::mutate(clstmemb = hv001.x == hv001.y) %>%
  dplyr::group_by(smpl1) %>%
  dplyr::summarise(
    n = n(),
    #clstprop = mean(clstmemb),
    F_wi = sum(malecotf)) %>% # \sum edge * malecotf (NB each edge is 1)
  dplyr::rename(name = smpl1)



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


node.edge.dens.weighted <- dplyr::left_join(node.edge.dens.weighted,
                                            mtdt, by = "name")

node.edge.dens.weighted %>%
  ggplot() +
  geom_point(aes(x = urban, y = F_wi))


node.edge.dens.weighted %>%
  dplyr::filter(F_wi > 49) %>%
  ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  geom_point(aes(x=longnum, y=latnum, color = F_wi), alpha = 0.5) +
  scale_color_viridis_c()

ggplot() +
  ggspatial::layer_spatial(data = urban, aes(fill = stat(band1))) +
  geom_sf(data = DRCprov, color = "#737373", fill = NA, size = 0.05) +
  scale_fill_distiller("Urbanicity", type = "div", palette = "RdYlBu")


summary(node.edge.dens.weighted$urban)


#..............................................................
# diff
#..............................................................
urban.pts <- mtdt %>%
  dplyr::select(c("name", "urban")) %>%
  dplyr::rename(smpl1 = name,
                urbanmean = urban)

ibD.urban <- dplyr::left_join(ibD, urban.pts, by = "smpl1")


colnames(urban.pts)[1] <- "smpl2"
ibD.urban <- dplyr::left_join(ibD.urban, urban.pts, by = "smpl2")

ibD.urban <- ibD.urban %>%
  dplyr::mutate(
    urbandiff = (urbanmean.x - urbanmean.y)^2
  )



ibD.urban.long <- expand_distance_matrix(ibD.urban) %>%
  dplyr::select(c("smpl1", "smpl2", "malecotf", "urbandiff"))
node.edge.dens.weighted <- ibD.urban.long %>%
  dplyr::mutate(F_wi = urbandiff*malecotf) %>%
  dplyr::group_by(smpl1) %>%
  dplyr::summarise(
    n = n(),
    F_wi = sum(F_wi)) %>% # \sum edge * malecotf (NB each edge is 1)
  dplyr::rename(name = smpl1)

node.edge.dens.weighted %>%
  ggplot() +
  geom_point(aes(x= , y = F_wi))





nick = ibD.urban %>%
  dplyr::filter(malecotf >= 0.5) %>%
  dplyr::filter(hv001.x != hv001.y)



