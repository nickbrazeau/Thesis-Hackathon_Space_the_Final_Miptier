#########################################################################
# Purpose: Pseudo-MLE to infer cluster level relatedness
#
# Author: Nicholas F. Brazeau
#
# Date: March 20 2020
#########################################################################

#' @title Expected Distribution
rij_expectation <- function(fi, fj, dij, m){
  ret <- ((fi + fj)/2) * exp(-dij*m)
  return(ret)
}


#' @title MSE Loss Function
#' @importFrom magrittr %>%
mse_loss_from_proposal <- function(yij, fij.prop, m.prop){
  # yij must be data.table
  # fij is a named vector, where names correspond to locations
  # set up key/dictionary approach, where rij is the "dict" and fij is the "key"

  #..............................................................
  # find proposal loss
  #..............................................................
  rij <- data.table::copy(yij) # have to be explicit to copy in place for DT (otherwise it acts as a pointer) https://stackoverflow.com/questions/10225098/understanding-exactly-when-a-data-table-is-a-reference-to-vs-a-copy-of-another
  rij[, `:=`(locat1 = as.character(locat1),
             locat2 = as.character(locat2)) ] # safety -- in case locations are names/numbers
  for(i in 1:length(fij.prop)) {
    rij[locat1 == names(fij.prop)[i], locat1 := fij.prop[i]]
    rij[locat2 == names(fij.prop)[i], locat2 := fij.prop[i]]
  }
  rij[, `:=`(locat1 = as.numeric(locat1),
             locat2 = as.numeric(locat2)) ] # fij are numerics, fix back

  # update expectation
  yij[, gendist_pred := rij_expectation(fi = rij[,locat1],
                                        fj = rij[,locat2],
                                        dij = rij[,geodist], m = m.prop) ]
  # calculate loss
  yij[, loss := (gendist - gendist_pred)^2 ]

  # estimate loss by location -- symmetrical
  ab.loss <- yij[, .(abloss = sum(loss)), by = locat1] %>%
    data.table::setnames(c("location", "abloss"))
  ba.loss <- yij[, .(baloss = sum(loss)), by = locat2] %>%
    data.table::setnames(c("location", "baloss"))
  ij.loss <- merge(ab.loss, ba.loss, by = c("location")) %>%
    .[, location_loss := ab.loss[, abloss] + ba.loss[,baloss]] %>%
    .[, c("location", "location_loss")]

  # return
  return(ij.loss)
}




#' @title Identify Relatedness in Continuous Space from Pairwise Comparisons
#'
#' @param gen.geo.dist dataframe; A dataframe with the following columns: Sample 1 Name; Sample 2 Name; Sample 1 Location; Sample 2 Location; Pairwise Genetic Distance; ; Pairwise Geographpic Distance
#' @param m_f_tuner numeric; weighted coin to see if M or F should be tuned
#' @param mexpshape numeric; gamma parameter in the exponential distribution for draws of M

relatedness_search <- function(gen.geo.dist, m_f_tuner = 0.5, mexpshape = 1, iters = 1e2){

  #..............................................................
  # Catches
  #..............................................................


  #..............................................................
  # Setup
  #..............................................................
  m.accept <- 0 # storage for moves
  fi.accept <- 0 # storage for moves
  loss.store <- rep(NA, length(iters))
  # for speed, let's go to data.table
  # for convenience, let's use column names
  gen.geo.dist <- data.table::as.data.table(gen.geo.dist)
  colnames(gen.geo.dist) <- c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist")

  # locations
  locations <- unique(unlist(gen.geo.dist[,c("locat1", "locat2")]))
  locations.update.counts <- matrix(NA, nrow = 1, ncol = length(locations))

  # standardize distance for stability
  gen.geo.dist[, geodist :=  geodist + 1e-7] # small offset
  distsd <- gen.geo.dist[, sd(geodist)]
  gen.geo.dist[, geodist :=  geodist /distsd]


  #..............................................................
  # Empirical Bayes Approach
  #..............................................................
  # inpsired by http://varianceexplained.org/r/empirical_bayes_baseball/
  # make prior based on whole dataset
  # beta is conjugate of binomial
  # success defined as any nonzero connection

  #..................
  # Global Estimate
  #..................
  # total pairwise counts -- remove self identifier chance, but this is the "count" of one slice of the square matrix
  totalpairwise <- length(unique(c(gen.geo.dist$smpl1, gen.geo.dist$smpl2)))
  totalpairwise <- totalpairwise^2 - totalpairwise

  # connections by location
  location1.connections <- gen.geo.dist[, c("smpl1", "locat1", "locat2", "gendist")]
  location2.connections <- gen.geo.dist[, c("smpl2", "locat1", "locat2", "gendist")]
  colnames(location2.connections) <-  c("smpl1", "locat2", "locat1", "gendist")
  priors <- rbind.data.frame(location1.connections, location2.connections) %>%
    .[, c("smpl1", "locat1", "gendist")]

  priors <- priors[, `:=`(success = sum(gendist != 0),
                          clustsmpls = length(unique(smpl1))),
                   by = locat1
                   ]
  priors <- unique(priors[, c("locat1", "success", "clustsmpls")])

  # sort by location name
  setkey(priors, locat1)

  # fit beta model (note, not doing a beta-binomial here, so not account for more weight on clusters with more obs... iid issue)
  mod.global <- fitdistrplus::fitdist(data =  priors$success/ (priors$clustsmpls * totalpairwise),
                                      distr = "beta",
                                      method = "mme", # use method of moments, which may be a bit more robust
                                      start = list(shape1 = 0.5, shape2 = 0.5))

  alpha.global <- mod.global$estimate[[1]]
  beta.global <- mod.global$estimate[[2]]


  #..................
  # Draw Fi from Prior
  # (Posterior) Distributions
  #..................
  # conjugate -- https://en.wikipedia.org/wiki/Conjugate_prior#Example
  priors[, `:=`(alpha1 = alpha.global + success,
                beta1  =  beta.global + totalpairwise - success)
         ]

  priors[, `:=`(mean = alpha1/(alpha1 + beta1),
                var  =  (alpha1*beta1)/( (alpha1 + beta1)^2 * (alpha1 + beta1 + 1) )
                )
         ]

  #..............................................................
  # INITIAL
  #..............................................................
  # Fij & M init prop
  fij.prop <- apply(priors, 1, function(x){return(rbeta(n = 1, shape1 = x["alpha1"], shape2 = x["beta1"]))})
  m.prop <- rexp(n = 1, rate = mexpshape)

  # store these for weights
  previous.weights <- matrix(1, nrow(priors), 1)
  previous.fis <- matrix(fij.prop, nrow(priors), 1)


  # observed
  yij <- gen.geo.dist[, c( "locat1", "locat2", "gendist", "geodist")]

  #...............
  # Calculate Loss from init
  #...............
  names(fij.prop) <- locations
  loss.curr <- mse_loss_from_proposal(yij, fij.prop, m.prop = m.prop)
  # store for out
  loss.store[1] <- sum(loss.curr[,location_loss])
  # store for rest
  fij.curr <- fij.prop
  m.curr <- m.prop

  #..............................................................
  # REST
  #..............................................................

  for (i in 2:iters){
  #..................
  # update M or update Fi's
  #..................

    if (rbinom(1,1, prob = m_f_tuner)) {
    #..................
    # update M
    #..................
      m.prop <- rexp(n = 1, rate = mexpshape)

      # Calculate loss and move
      loss.prop <- mse_loss_from_proposal(yij, fij.curr, m.prop = m.prop)
      if (sum(loss.prop$location_loss) < sum(loss.curr$location_loss)) {
        # save out new items
        m.accept <- m.accept + 1
        m.curr <- m.prop
        loss.curr <- loss.prop
      }
    } else {
    #..................
    # update Fis
    #..................
      fij.prop <- apply(priors, 1, function(x){return(rbeta(n = 1, shape1 = x["alpha1"], shape2 = x["beta1"]))})
      names(fij.prop) <- priors$locat1

      # Calculate loss and move
      loss.prop <- mse_loss_from_proposal(yij, fij.prop, m.prop = m.curr)
      if (sum(loss.prop$location_loss) < sum(loss.curr$location_loss)) {
        fi.accept <- fi.accept + 1

        #..................
        # update Priors based on loss
        #..................
        previous.fis <- cbind(previous.fis, fij.prop)
        newwi <- 1/(loss.prop$location_loss/sum(loss.prop$location_loss)) # inverse of how far off, so bigger loss is smaller number
        newwi <- newwi/sum(newwi) # make it a weight
        previous.weights <- cbind(previous.weights, newwi)
        previous.weights[,1] <- previous.weights[,1] * 1/(ncol(previous.weights)-1) # dampening factor once N is large
        newmeans <- rowSums(previous.fis * previous.weights)/rowSums(previous.weights)

        # the mean and sample size parameters are related to the shape parameters α and β via α = μν, β = (1 − μ)ν
        # this isn't haldane's prior though, so this isn't proper bayesian
        priors[, `:=`(alpha1 = newmeans * (alpha1+beta1),
                      beta1 = (1-newmeans) * (alpha1+beta1))]


        # save out new items
        fi.accept <- fi.accept + 1
        fij.curr <- fij.prop
        loss.curr <- loss.prop

      }
    } # end loop for fi

    # store loss
    loss.store[i] <- sum(loss.curr[,location_loss])
  }

  #..............................................................
  # return
  #..............................................................
  posteriors <- priors[,c("locat1", "alpha1", "beta1")]
  posteriors[, `:=`(mean = alpha1/(alpha1 + beta1),
                    var = (alpha1*beta1)/( ((alpha1 + beta1)^2) * (alpha1 + beta1 + 1) )
                    )]

  ret <- list(
    fi_betas = as.data.frame(posteriors), # drop data table class
    migration = m.curr,
    loss.store = loss.store,
    m.accept = m.accept,
    fi.accept = fi.accept
    )
  return(ret)

}
