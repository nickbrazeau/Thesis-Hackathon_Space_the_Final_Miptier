#............................................................
# should i look to see if coercing the continuous spatial process
# to discrete causes issues?
#........................................................
#and nick do a pcoa...

## .................................................................................
## Purpose: Look at other markers of diversity and gene flow
##
## Notes:
## .................................................................................
library(tidyverse)
library(adegenet)
library(hierfstat)
library(vcfR)
source("R/themes.R")

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
keda <- data.frame(k = seq(2, 200, by = 1))
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
  geom_vline(xintercept = 75, color = "red", linetype = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "blue", linetype = 2, alpha = 0.8) +
  theme_minimal() +
  ylab("Total Within-Cluster Sum of Squares") +
  xlab("K")

kesteda.plotObj

#......................
# define population partitions
#......................
k <- kmeans(coords[,2:3], 75)
demes <- cbind.data.frame(coords, k = k$cluster)
deme_mtdt <- mtdt %>%
  dplyr::select(c("name", "hv001")) %>%
  dplyr::left_join(., demes, by = "hv001") %>%
  dplyr::mutate(k = factor(k))

#......................
# viz parts
#......................
deme_mtdt %>%
  dplyr::group_by(k) %>%
  dplyr::summarise(
    cnt = dplyr::n(),
    longnum = mean(longnum),
    latnum = mean(latnum)
  ) %>%
  ggplot() +
  geom_point(aes(x = longnum, y = latnum,
                 size = cnt, color = cnt),
             alpha = 0.7) +
  scale_color_viridis_c() +
  ggtitle("By K")


deme_mtdt %>%
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    k = unique(k),
    longnum = mean(longnum),
    latnum = mean(latnum)
  ) %>%
  ggplot() +
  geom_point(aes(x = longnum, y = latnum,
                 color = k),
             alpha = 0.7) +
  ggtitle("By hv001") +
  theme(legend.position = "none")


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

#......................
# get pairwise Fst w/ each deme as outgroup
#......................
get_outgroup_fst <- function(x, gnind_k){
  change <- which(gnind_k@pop != x)
  # refactor
  gnind_k@pop <- factor(gnind_k@pop, levels = c(levels(gnind_k@pop), "out"))
  gnind_k@pop[change] <- "out"
  dat <- hierfstat::genind2hierfstat(gnind_k)
  wcFst <- hierfstat::pairwise.WCfst(dat, diploid = T)
  return(wcFst[lower.tri(wcFst)])
}
# outgroup fst values
wcFst_outgroup <- tibble::tibble(kvals = unique(deme_mtdt$k))
wcFst_outgroup <- wcFst_outgroup %>%
  dplyr::mutate(fstout = map_dbl(kvals, get_outgroup_fst, gnind_k = gnind))

# store
Fstout_deme_mtdt <- mtdt %>%
  dplyr::select(c("name", "hv001")) %>%
  dplyr::left_join(., demes, by = "hv001") %>%
  dplyr::mutate(kvals = factor(k)) %>%
  dplyr::left_join(., wcFst_outgroup, by = "kvals") %>%
  dplyr::mutate(fstout = ifelse(fstout < 0, 0, fstout))

# save out
dir.create("results/gendiversity/", recursive = T)
saveRDS(Fstout_deme_mtdt, "results/gendiversity/fst_outgroup_ests.RDS")

# viz out
Fstout_deme_mtdt %>%
  ggplot() +
  geom_point(aes(x = longnum, y = latnum,
                 color = fstout),
             alpha = 0.7) +
  scale_color_viridis_c()


#............................................................
# Population Nucleotide Diversity
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
# Visualizing Wave Forms of WSAFs
#   this is more or less what is under the hood
#   of a gaussian mixture model
#...........................................................
