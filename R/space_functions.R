#' @title Krigging by fields package
#' @description function is not generalizable and is meant for this project use case
get_krigged <- function(discentret, locats, range, res) {
  # check
  if(!all(colnames(locats) %in% c("longnum", "latnum", "deme"))){
    stop("locats object not right")
  }

  # make discent return object
  ret <- tibble::tibble(disc = discentret$Final_Fis,
                        key = discentret$deme_key$key,
                        deme = discentret$deme_key$Deme) %>%
    dplyr::left_join(., locats, by = "deme")

  #......................
  # krig function
  #   https://github.com/NCAR/fields/blob/master/fieldsVignette.pdf
  #......................

  # find an appropriate theta
  mle_fits <- fields::spatialProcess(cbind(ret$longnum, ret$latnum),
                                     y = ret$disc,
                                     cov.args = list(Covariance = "Matern",
                                                     smoothness = 0.5)) # Note that the exponential is also part of the Matern family with smoothness set to .5

  # fit model
  krig_fit <- fields::Krig(x = cbind(ret$longnum, ret$latnum),
                           Y = ret$disc,
                           theta = mle_fits$theta.MLE)

  # predict at new places
  krig_preds <- predictSurface(object = krig_fit,
                               nx = res, ny = res)
  krig_ses <- predictSurfaceSE(object = krig_fit,
                                 nx = res, ny = res)


  #......................
  # out
  #......................
  # bring together
  preds <- cbind.data.frame(expand.grid(krig_preds$x, krig_preds$y),
                          as.vector(krig_preds$z))
  colnames(preds) <- c("longnum", "latnum", "pred_disc")

  ses <- cbind.data.frame(expand.grid(krig_ses$x, krig_ses$y),
                            as.vector(krig_ses$z))
  colnames(ses) <- c("longnum", "latnum", "pred_disc_ses")
  # out
  out <- dplyr::left_join(preds, ses, by = c("longnum", "latnum"))

  return(tibble::as_tibble(out))


}


#' #' @title Krigging by Verity
#' #' @description function is not generalizable and is meant for this project use case
#'
#' get_krigged <- function(discentret, locats, new_long_lat) {
#'   # check
#'   if(!all(colnames(locats) %in% c("longnum", "latnum", "deme"))){
#'     stop("locats object not right")
#'   }
#'   if(c(!all(colnames(new_long_lat) %in% c("longnum", "latnum"))) &
#'      !is.matrix(new_long_lat)){
#'     stop("new_long_lat object not right")
#'   }
#'
#'   # make discent return object
#'   ret <- tibble::tibble(disc = discentret$Final_Fis,
#'                         key = discentret$deme_key$key,
#'                         deme = discentret$deme_key$Deme) %>%
#'     dplyr::left_join(., locats, by = "deme")
#'
#'   #......................
#'   # do krig w/ verity neg like and gp approach
#'   #......................
#'   # calculate -log(likelihood) given GP parameters
#'   negLogLike <- function(params, distMat_squared, z) {
#'     s2_f <- params[1]
#'     s2_n <- params[2]
#'     l <- params[3]
#'     K <- s2_f*exp(-distMat_squared/(2*l))
#'     K <- K + diag(s2_n,length(z))
#'     K_det <- determinant(K,logarithm=TRUE)$modulus[1]
#'     L <- chol(K)
#'     logLike <- -0.5*t(z)%*%backsolve(L, forwardsolve(t(L),z)) - 0.5*log(K_det)
#'     return(-logLike)
#'   }
#'
#'   # estimate GP parameters by maximum likelihood
#'   distMat_squared <- as.matrix(dist(cbind(ret$longnum, ret$latnum)))^2
#'   opt <- optim(par = c(1,1,1),
#'                distMat_squared = distMat_squared,
#'                z = ret$disc,
#'                fn = negLogLike,
#'                method="L-BFGS-B")
#'   s2_f <- opt$par[1]
#'   s2_n <- opt$par[2]
#'   l <- opt$par[3]
#'
#'   # estimate mean of GP
#'   K <- s2_f*exp(-distMat_squared/(2*l))
#'   K <- K + diag(s2_n, length(ret$longnum))
#'   K_inv <- solve(K)
#'   distMat_new <- fields::rdist(new_long_lat, cbind(ret$longnum, ret$latnum))
#'   Ks <- s2_f*exp(-distMat_new^2/(2*l))
#'   z_new <- Ks%*%K_inv%*%ret$disc
#'
#'   #......................
#'   # out
#'   #......................
#'   out <-cbind.data.frame(new_long_lat, z_new)
#'   colnames(out) <- c("longnum", "latnum", "pred_disc")
#'
#'   return(tibble::as_tibble(out))
#' }


#' @title Basic Point Plotter
#' @description super simple function is not generalizable and is meant for this project use case
raster_plotter <- function(pred_discentret) {
  if(!all(colnames(pred_discentret) %in% c("longnum", "latnum", "pred_disc"))){
    stop("pred_discentret object not right")
  }

  # make raster
  rstrdat <- raster::rasterFromXYZ(xyz = pred_discentret)
  # make plot
  plotObj <- ggplot() +
    ggspatial::layer_spatial(data = rstrdat, aes(fill = stat(band1))) +
    geom_contour(data = pred_discentret,
                 aes(x = longnum, y = latnum, z = pred_disc), color = "black")
  # out and let user do graph things (eg themes)

  return(plotObj)
}



#' @title Basic Point Plotter
#' @description plot points for disc
point_plotter <- function(discentret, locats) {
  if(!all(colnames(locats) %in% c("longnum", "latnum", "deme"))){
    stop("locats object not right")
  }

  # make discent return object
  ret <- tibble::tibble(disc = discentret$Final_Fis,
                        key = discentret$deme_key$key,
                        deme = discentret$deme_key$Deme)

  # put back in geo locations
  plotObj <- dplyr::left_join(ret, locats, by = "deme") %>%
    ggplot() +
    geom_point(aes(x = longnum, y = latnum, color = disc),
               size = 1.2, alpha = 0.9)
  # out and let user do graph things (eg themes)
  return(plotObj)
}


#' @title DISC Estimates & Realized IBD
#' @description Combine disc estimates and realized ibd plots
plot_realIBD_discent <- function(gengeodat, discentret, locats) {
  if(!all(colnames(locats) %in% c("longnum", "latnum", "deme"))){
    stop("locats object not right")
  }

  #......................
  # DISC estimates
  #......................
  # make discent return object
  ret <- tibble::tibble(disc = discentret$Final_Fis,
                        key = discentret$deme_key$key,
                        deme = discentret$deme_key$Deme)

  # put back in geo locations
  discplotdat <- dplyr::left_join(ret, locats, by = "deme")

  #......................
  # realized ibd
  #......................
  # quick manips
  locats1 <- locats %>%
    dplyr::rename(locat1 = deme)
  locats2 <- locats %>%
    dplyr::rename(locat2 = deme)
  # bring together realized
  realIBD_plotdat <- gengeodat %>%
    dplyr::left_join(., y = locats1, by = "locat1") %>%
    dplyr::left_join(., y = locats2, by = "locat2") %>%
    dplyr::filter(gendist > 0 )

  #......................
  # bring together and out
  #......................
  ggplot() +
    geom_point(data = discplotdat, aes(x = longnum, y = latnum, fill = disc),
               size = 4, alpha = 1, shape = 22, color = "#000000") +
    geom_segment(data = realIBD_plotdat,
                 aes(x = longnum.x,
                     y = latnum.x,
                     xend = longnum.y,
                     yend = latnum.y,
                     color = gendist),
                 alpha = 0.5) +
    scale_color_viridis_c("Sim. IBD", option = "viridis") +
    scale_fill_viridis_c("Pred. DISC", option = "cividis") +
    theme_void()

}
