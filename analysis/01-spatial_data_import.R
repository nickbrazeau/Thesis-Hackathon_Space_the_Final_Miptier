#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import spatial data
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)

gdrive <- tcltk::tk_choose.dir()


#.............
# GADM Downloadds
#.............
#spatial from GADM -- these are polygon files, doing this is legacy as raster does this nicely...but already fixed naming issue here
DRCprov <- sf::st_as_sf( raster::getData(name = "GADM", country = "CD", level = 1, path = paste0(gdrive, "/data/map_bases/gadm/")) )
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

saveRDS(DRCprov, file = paste0(gdrive, "/data/map_bases/cd2013_drcprov.rds"))








#.............
# Accessibility from Weiss
#.............
# see explanation in the DHS GC manual
# NOTE, this is from 2015
summary(dt$travel_times_2015)
hist(dt$travel_times_2015)
dt <- dt %>%
  dplyr::mutate(travel_times_2015_scale = scale(log(travel_times_2015 + tol), center = T, scale = T))

summary(dt$travel_times_2015_scale)
hist(dt$travel_times_2015_scale) # many, many 0s
hist(dt$travel_times_2015_scale) # standardization doesn't look as good, but should capture urban v. rural well
summary(dt$travel_times_2015_scale); sd(dt$travel_times_2015_scale)


