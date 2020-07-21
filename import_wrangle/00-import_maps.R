#----------------------------------------------------------------------------------------------------
# Purpose of this script is to download map features that will be need for later plotting
# Will use this to mostly make "pretty" maps
# http://www.francescobailo.net/2018/08/how-to-quickly-enrich-a-map-with-natural-and-anthropic-details/
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(sf)
library(raster)
library(ggspatial)
source("R/basics.R")
dir.create("data/map_bases/gadm/", recursive = T)

#---------------------------------------------------------------------------------
# pull down DRC maps from GADM
#---------------------------------------------------------------------------------
#spatial from GADM -- these are polygon files, doing this is legacy as raster does this nicely...but already fixed naming issue here
if(!dir.exists(paste0(getwd(), "/data/map_bases/gadm/"))){dir.create(paste0(getwd(), "/data/map_bases/gadm/"), recursive = T)}
DRCprov <- sf::st_as_sf( raster::getData(name = "GADM", country = "CD", level = 1, path = "data/map_bases/gadm/") )
colnames(DRCprov) <- tolower(colnames(DRCprov))
colnames(DRCprov)[4] <- "adm1name" # to match the DHS province names
# need to strip accent marks also to match the DHS province names
# https://stackoverflow.com/questions/20495598/replace-accented-characters-in-r-with-non-accented-counterpart-utf-8-encoding
# thanks to @Thomas for this great trick
unwanted_array = list(   'Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 'Ç'='C', 'È'='E', 'É'='E',
                         'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U',
                         'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 'ß'='Ss', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c',
                         'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
                         'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y' )

DRCprov$adm1name <- chartr(paste(names(unwanted_array), collapse=''),
                           paste(unwanted_array, collapse=''),
                           DRCprov$adm1name)

# match dhs and gadm
DRCprov$adm1name[DRCprov$adm1name == "Kongo-Central"] <- "Kongo Central"
DRCprov$adm1name[DRCprov$adm1name == "Tanganyika"] <- "Tanganyka" # note, DHS has the typo here

# overwrites
saveRDS(DRCprov, "data/map_bases/gadm/gadm36_COD_1_sp.rds")


#..............................
# Pull down the border cntrs w/ raster
#..............................
brdrcnt <- lapply(c("UGA", "SSD", "CAF", "COG", "AGO", "ZMB", "TZA", "RWA", "BDI", "GAB", "CMR", "GNQ"),
                  function(x){
                    ret <- raster::getData(name = "GADM", country = x, level = 0, path = "data/map_bases/gadm/")
                    ret <- sf::st_as_sf(ret)
                    return(ret)

                  })

drcadm0 <- raster::getData(name = "GADM", country = "COD", level = 0, path = "data/map_bases/gadm/")

#..............................
# Pull down ocean
#..............................
# https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip
oceans <- sf::st_read("data/map_bases/ne_10m_ocean/ne_10m_ocean.shp")

#---------------------------------------------------------------------------------
# write out lists of map bases for later plotting
#---------------------------------------------------------------------------------

prettybasemap_nodrc <- list(
  geom_sf(data = brdrcnt[[1]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  # geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']),
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']),
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",
                                    pad_y = unit(1.25, "cm")),
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite  theme
)



prettybasemap_nodrc_nonorth <- list(
  geom_sf(data = brdrcnt[[1]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  # geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']),
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']),
           datum = NA),
  #   ggspatial::annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(1.25, "cm")),
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite  theme
)

prettybasemap_nodrc_dark <- list(
  geom_sf(data = brdrcnt[[1]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#525252", color = "#737373",lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']),
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']),
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",
                                    pad_y = unit(1.25, "cm")),
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite  theme
)

prettybasemap_nodrc_nonorth_dark <- list(
  geom_sf(data = brdrcnt[[1]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#525252", color = "#737373",lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#525252", color = "#737373", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']),
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']),
           datum = NA),
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite vivid theme
)





#----------------------------------------------------------------------------------------------------
# Save Objects & Write out
#----------------------------------------------------------------------------------------------------
save(smpl_bckgrnd,
     prettybasemap_nodrc,
     prettybasemap_nodrc_nonorth,
     prettybasemap_nodrc_dark,
     prettybasemap_nodrc_nonorth_dark,
     file = "data/map_bases/space_mips_maps_bases.rda")
saveRDS(DRCprov, file = "data/map_bases/spacemips_DRCprov.rds")





