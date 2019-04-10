#----------------------------------------------------------------------------------------------------
# Read Raster Temp and Precip
#----------------------------------------------------------------------------------------------------
readRasterBB <- function(rstfile, bb = bb){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = sf::as_Spatial(bb))
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m

  return(ret)

}

extractDHSpts_from_raster <- function(hv001, raster,
                                      geometry, urban_rura){
  # note raster is expected to have units of m

  geometry <- sf::as_Spatial(geometry)
  buffer = ifelse(urban_rura == "R", 10, 2)
  ret <- raster::extract(x = raster[[1]],
                         y =  geometry[[1]],
                         buffer = buffer,
                         fun = mean,
                         sp = F)
  return(ret)
}


