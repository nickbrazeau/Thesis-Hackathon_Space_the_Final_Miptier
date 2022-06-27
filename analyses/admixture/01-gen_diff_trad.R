## .................................................................................
## Purpose: Look at other markers of diversity and gene
## flow using traditional approaches
##
## Notes:
## .................................................................................
library(tidyverse)
library(adegenet)
library(hierfstat)
library(vcfR)
library(ape)
source("R/themes.R")
DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))

#............................................................
# tidy data
#...........................................................
#......................
# read in data
#......................
# meta
mtdt <- readRDS("data/derived_data/sample_metadata.rds")
# genetic
vcf <- vcfR::read.vcfR("data/derived_data/bigbarcode_genetic_data/mipbivcfR.DRC.vcf.gz")
gnind <- vcfR::vcfR2genind(vcf)

# make sure we get names in right order
adgnm <- tibble::tibble(name = indNames(gnind))

#............................................................
# K-means for populations
#...........................................................
coords <- mtdt %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>%
  dplyr::select(c("hv001", "longnum", "latnum")) %>%
  dplyr::filter(!duplicated(.))

#.............
# K-means Clustering of DRC
#.............
keda <- data.frame(k = seq(2, 75, by = 1))
keda$kmeans <- map(keda$k, function(k){return(kmeans(x = coords[,2:3], centers = k))})
keda$wss <- map(keda$kmeans, "withinss")
keda$totalwss <- map(keda$wss, function(x){return(sum(x))})

keda.df <- keda %>%
  dplyr::select(c("k", "totalwss")) %>%
  tidyr::unnest(cols = c("totalwss"))

kesteda.plotObj <- keda.df %>%
  tibble::as_tibble(.) %>%
  ggplot() +
  geom_line(aes(x=k, y=totalwss)) +
  geom_point(aes(x=k, y=totalwss)) +
  geom_vline(xintercept = 15, color = "red", linetype = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "blue", linetype = 2, alpha = 0.8) +
  theme_minimal() +
  ylab("Total Within-Cluster Sum of Squares") +
  xlab("K")

kesteda.plotObj



#......................
# define population partitions
#......................
k <- kmeans(coords[,2:3], 15)
demes <- cbind.data.frame(coords, k = k$cluster)
deme_mtdt <- mtdt %>%
  dplyr::select(c("name", "hv001")) %>%
  dplyr::left_join(., demes, by = "hv001") %>%
  dplyr::mutate(k = factor(k))

# quick look at deme order
deme_mtdt %>%
  ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(aes(x = longnum, y = latnum, color = factor(k)), show.legend = F) +
  scale_color_manual("", values = c(brewer.pal(n = 8, name = "Dark2"),
                                    "#3C3B6E", "#B22234")) +
  theme_void()

#......................
# combine in adgenet
#......................
# bring in adegenet names
adgnm <- dplyr::left_join(adgnm, deme_mtdt, by = "name")
# set adegent population
gnind@pop <- adgnm$k

#............................................................
# diversity metrics
#...........................................................
LociHe <- summary(gnind)
# loci diversity
plot(LociHe$Hobs, LociHe$Hexp, xlab="Hobs", ylab="Hexp",
     main="Expected heterozygosity as a function of observed heterozygosity per locus")
div <- hierfstat::genind2hierfstat(gnind) # Create hierfstat object
basestat <- hierfstat::basic.stats(div, diploid = T)

# deme size
wcFst <- hierfstat::pairwise.WCfst(div, diploid = T)

# make plot
wcFstdf <- broom::tidy(as.dist(wcFst)) %>%
  dplyr::mutate(item1 = forcats::fct_rev(forcats::fct_reorder(.f = item1, .x = item1, .fun = length)),
                item2 = forcats::fct_rev(forcats::fct_reorder(.f = item2, .x = item2, .fun = length)))
# heatplot
wcFstdf %>%
  ggplot() +
  geom_tile(aes(x = item1, y = item2, fill = distance)) +
  scale_fill_viridis_c("Fst")

#......................
# make panel A
#......................
mn <- wcFstdf %>%
  ggplot() +
  geom_tile(aes(x = item1, y = item2, fill = distance)) +
  scale_fill_viridis_c("Fst") +
  plot_theme +
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.25, vjust = 0.5))

inset <- deme_mtdt %>%
  ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(aes(x = longnum, y = latnum, color = factor(k)),
             size = 0.5) +
  scale_color_manual("Pop.", values = c(brewer.pal(n = 8, name = "Dark2"),
                                        "#C4B274", "#B22234",
                                        "#0D374F", "#48736C",
                                        "#AFC3AE", "#703837", "#966036")) +
  theme_void() +
  theme(legend.title = element_text(face = "bold", hjust = 0.5, vjust = 0.5),
        legend.text = element_text(face = "bold")) +
  guides(color = guide_legend(ncol = 2))



# bring together
panelA <- cowplot::ggdraw() +
  cowplot::draw_plot(mn,
                     x = 0, y = 0, width = 1, height = 0.95, scale = 1) +
  cowplot::draw_plot(inset, x = 0.5, y = 0.45,
                     width = 0.5, height = 0.5)



#............................................................
# Population Nucleotide Diversity
#   not much signal here
#...........................................................
popdiv <- adegenet::Hs(gnind)
summary(popdiv)
table(round(popdiv, digits = 2))

demes %>%
  dplyr::mutate(k = as.character(k)) %>%
  dplyr::left_join(., tibble::tibble(k = names(popdiv),
                                     div = unname(popdiv))) %>%
  ggplot() +
  geom_point(aes(x = longnum, y = latnum,
                 color = div),
             alpha = 0.7) +
  scale_color_distiller("Diversity", type = "div", palette = "RdYlBu",
                        na.value = NA) +
  map_theme +
  theme(plot.background = element_rect(fill = "#000000"))


#............................................................
# PCOA
#...........................................................
#......................
# read in distance calculation from Amato et al (PMC4786412)
#......................
gendiff <- readr::read_tsv("analyses/admixture/mipbivcfR_DRC_dab.tab.txt",
                           col_names = F) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "numsites", "dab", "dab_wi"))
# make diagonal for matrix
gendiag <- tibble::tibble(smpl1 = unique(c(gendiff$smpl1, gendiff$smpl2)),
                          smpl2 = unique(c(gendiff$smpl1, gendiff$smpl2)),
                          numsites = NA,
                          dab = 1,
                          dab_wi = 1)
# add diagonal and order
gendiff <- dplyr::bind_rows(gendiff, gendiag) %>%
  dplyr::mutate(smplfct1 = forcats::fct_rev(forcats::fct_reorder(.f = smpl1, .x = smpl1, .fun = length)),
                smplfct2 = forcats::fct_rev(forcats::fct_reorder(.f = smpl2, .x = smpl2, .fun = length)))

gendiffmatdf <- gendiff %>%
  dplyr::select(c("smplfct1", "smplfct2", "dab_wi")) %>%
  tidyr::pivot_wider(data = ., names_from = smplfct2, values_from = dab_wi) %>%
  dplyr::select(c("smplfct1", "1023W9S2B", dplyr::everything()))
# conver to matrix
gendiffmat <- as.matrix(gendiffmatdf[,2:ncol(gendiffmatdf)])
rownames(gendiffmat) <- gendiffmatdf$smplfct1
# sanity check
sum(diag(gendiffmat) != 1)
# get lower triangle
gendiffmat[lower.tri(gendiffmat)]  <- t(gendiffmat)[lower.tri(gendiffmat)]

#......................
# run ape pcoa
#......................
retpcoa <- ape::pcoa(gendiffmat)
# get eigenvales
eigendf <- tibble::as_tibble(retpcoa$vectors, .name_repair = "minimal") %>%
  dplyr::mutate(name = rownames(retpcoa$vector)) %>%
  dplyr::select(c("name", dplyr::everything())) %>%
  magrittr::set_colnames(c("name", "pe1", "pe2", "pe3", "pe4"))

# bring in metadata
dplyr::left_join(eigendf, mtdt) %>%
  dplyr::select(c("name", "longnum", "latnum", dplyr::starts_with("pe"))) %>%
  tidyr::pivot_longer(data = ., cols = -c("name", "longnum", "latnum"),
                      names_to = "eigenaxis", values_to = "poseigen") %>%
  ggplot() +
  geom_point(aes(x = longnum, y = latnum, color = poseigen)) +
  facet_wrap(~eigenaxis, scales = "free") +
  scale_color_viridis_c()


dplyr::left_join(eigendf, mtdt) %>%
  ggplot() +
  geom_point(aes(x = longnum, y = latnum, color = pe1)) +
  scale_color_viridis_c()


#......................
# get radial distance
#......................
DRC <- sf::st_union(DRCprov)
# center
drc_center <- sf::st_centroid(DRC)
drc_center <- sf::st_coordinates(drc_center)
# simple euclidean distance
eigendfmtdt <- dplyr::left_join(eigendf, mtdt)
eigendfmtdt_sf <- sf::st_as_sf(eigendfmtdt, coords = c("longnum", "latnum"),
                               crs = "+init=epsg:4326")

eigendfmtdt <- eigendfmtdt %>%
  dplyr::mutate(raddist = as.numeric(sf::st_distance(eigendfmtdt_sf, drc_center))/1e3)

# sanity check
eigendfmtdt %>%
  ggplot() +
  geom_point(aes(x = longnum, y = latnum, color = raddist))


#......................
# use bob's scatterplot 3d
#......................
panelB <- bobfunctions2::gg3d_scatterplot(x = eigendf$pe1,
                                          y = eigendf$pe2,
                                          z = eigendf$pe3,
                                          colour = eigendfmtdt$raddist,
                                          size = 1.2,
                                          alpha = 0.8,
                                          theta = 135,
                                          phi = 15,
                                          d = 2,
                                          axis_on = T,
                                          z_type = 1,
                                          tick_length = 0.005,
                                          x_lim = c(-0.1, 0.03),
                                          y_lim = c(-0.1, 0.03),
                                          z_lim = c(-0.05, 0.03),
                                          axis_lab_size = 0,
                                          axis_lab_dist = 0) +
  scale_color_viridis_c("Radial Distance \nfrom DRC Centroid",
                        option = "cividis") +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold", angle = 45))


#............................................................
# bring it home for final figure
#...........................................................
jpeg("results/figures/investigate_admix_trad.jpg",
     width = 11, height = 6, units = "in", res = 500)
cowplot::plot_grid(panelA, panelB,
                   ncol = 2, nrow = 1,
                   labels = c("(A)", "(B)"),
                   rel_heights = c(0.5, 1))
graphics.off()
