library(ggplot2)

#----------------------------------------------------------------------------------------------------
# Basic epi
#----------------------------------------------------------------------------------------------------

logit <- function(x, tol=1e-4){
  return( log(((x+tol)/(1-x+tol))) )
}

#----------------------------------------------------------------------------------------------------
# Make final Survey Object to account for DHS Survey Weights
#----------------------------------------------------------------------------------------------------
makecd2013survey <- function(survey = dt){
  options(survey.lonely.psu="certainty")
  dtsrvy <- survey %>% srvyr::as_survey_design(ids = hv001,
                                               strata = hv023,
                                               weights = hv005_wi)
  return(dtsrvy)
}

#----------------------------------------------------------------------------------------------------
# Change return from scale
#----------------------------------------------------------------------------------------------------
my.scale <- function(x, ...){
  ret <- scale(x, ...)
  if(ncol(ret) == 1){
    return(as.numeric(ret)) # coerce matrix of 1 col to numeric vector
  } else{
    stop("Matrix has more than one column")
  }
}


