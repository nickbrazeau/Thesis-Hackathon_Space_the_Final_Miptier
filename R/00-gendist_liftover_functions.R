


## These  two functions need to be better/slicker
## What seems sensible is to cache the function call with the
## returned object so that when we call append_genetic_distances it checks
## to see if we have already calculated the distance and thus does not need to
## But I run out of time to make it fancy like this.


append_genetic_distances <- function(mipanalyzer_object,
                                     dist_func,
                                     ...) {
  # Checks
  # ----------------------------------------------------------------------------
  MIPanalyzer:::assert_custom_class(x, "mipanalyzer_biallelic")

  # Setting up variables for where to store results
  # ----------------------------------------------------------------------------

  # default distances (should perhaps incorporate the arguments so they are not
  # overwritten with different argument sets)
  distance_functions <- c(
    "IBS_major_allele" = "MIPanalyzer::get_IBS_distance",
    "WSAF_comparison" = "MIPanalyzer::get_IB_mixture",
    "DAB_genomic" = "MIPanalyzer::get_genomic_distance",
    "MLE_IBD" = "MIPanalyzer::inbreeding_mle"
  )

  # sort out which distance function we are using
  nm <- grep(names(dist_func), names(distance_functions))
  dist_func <- dist_func[[1]]

  # check we have that function
  if (length(nm) == 0) {
    stop("No matching distance function")
  } else {
    nm <- names(distance_functions[nm])
  }

  # create the gen_distances list if needed
  if (is.null(mipanalyzer_object$gen_distances)) {
    mipanalyzer_object$gen_distances <- list()
  }

  # Calculate distances
  # ----------------------------------------------------------------------------

  # calculate and store the results
  if (nm != "MLE_IBD") {
    mipanalyzer_object$gen_distances[[nm]] <-
      dist_func(mipanalyzer_object, ...)
  } else {
    mipanalyzer_object$gen_distances[[nm]] <-
      dist_func(mipanalyzer_object, ...)$mle
  }

  # attribute with function call (would be good to have the actual call)
  attr(mipanalyzer_object$gen_distances[[nm]], "function_call") <-
    distance_functions[nm]

  return(mipanalyzer_object)

}


add_genetic_distances <- function(mipanalyzer_object,
                                  args = NULL) {
  # Checks
  # ----------------------------------------------------------------------------
  MIPanalyzer:::assert_custom_class(x, "mipanalyzer_biallelic")

  # Setting up variables for where to store results
  # ----------------------------------------------------------------------------

  # defaults
  distance_functions <- c(
    "IBS_major_allele" = MIPanalyzer::get_IBS_distance,
    "WSAF_comparison" = MIPanalyzer::get_IB_mixture,
    "DAB_genomic" = MIPanalyzer::get_genomic_distance,
    "MLE_IBD" = MIPanalyzer::inbreeding_mle
  )

  # default function args
  if (is.null(args)) {
    args <- list(
      "IBS_major_allele" = list(ignore_het = FALSE),
      "WSAF_comparison" = list(tol = 0.05),
      "DAB_genomic" = list(cutoff = 0.1),
      "MLE_IBD" = list(f = seq(0.01, 0.99, 0.01), ignore_het = FALSE)
    )
  }

  # Calculate distances
  # ----------------------------------------------------------------------------

  # loop through distance functions and apply
  for (i in seq_len(length(args))) {
    mipanalyzer_object <-
      append_genetic_distances(mipanalyzer_object = mipanalyzer_object,
                               dist_func = distance_functions[i],
                               unlist(args[i]))
  }

  return(mipanalyzer_object)

}

add_admin_summaries <- function(mipanalyzer_object) {

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results <- list()
  }

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results$admin_gen_dist <- list()
  }

  distances <- names(mipanalyzer_object$gen_distances)

  for (i in distances) {
    mipanalyzer_object$results$admin_gen_dist[[i]] <-
      getadmin_gendist_summary(
        mipanalyzerobject_samples = mipanalyzer_object$samples,
        mipanalyzerobject_gendistmat = mipanalyzer_object$gen_distances[i],
        type = "mean" # TODO fix matcharg
      )
  }

  return(mipanalyzer_object)

}


add_cluster_summaries <- function(mipanalyzer_object) {

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results <- list()
  }

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results$cluster_gen_dist <- list()
  }

  distances <- names(mipanalyzer_object$gen_distances)

  for (i in distances) {
    mipanalyzer_object$results$cluster_gen_dist[[i]] <-
      get_cluster_summaries(
        mipanalyzerobject_samples = mipanalyzer_object$samples,
        mipanalyzerobject_gendistmat = mipanalyzer_object$gen_distances[i],
        type = "mean" # TODO fix matcharg
      )
  }

  return(mipanalyzer_object)

}

add_raw_summaries <- function(mipanalyzer_object) {

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results <- list()
  }

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results$raw_gen_dist <- list()
  }

  distances <- names(mipanalyzer_object$gen_distances)

  for (i in distances) {
    mipanalyzer_object$results$raw_gen_dist[[i]] <-
      get_raw_summaries(
        mipanalyzerobject_samples = mipanalyzer_object$samples,
        mipanalyzerobject_gendistmat = mipanalyzer_object$gen_distances[i],
        type = "mean" # TODO fix matcharg
      )
  }

  return(mipanalyzer_object)

}


# The following summary functions should be tidied up to better drop in
# which other covariates you want to group in etc in the key creation etc
# as well as the type argument to generecise the summary stat creation
# -----------------------------------------------------------------------

# Rationale in how these are summarised with respect to adding coavraiates later is:

# For admins we will just left bind admin summaries saved at paste0(gdrive, "/data/derived_data/cd2013_kids_dhs_admin1_recode.rds") (admin summaries have _ in this rds)
# For cluster we will just left bind cluster summaries saved at paste0(gdrive, "/data/derived_data/cd2013_kids_dhs_recode.rds") (cluster summaries have _ in this rds)
# For individual we will just left bind raw individual data saved at paste0(gdrive, "/data/derived_data/cd2013_kids_dhs_recode.rds") (i.e. vars without _ in this rds)

# summarise gentic distance by admin
get_admin_gendist_summary <- function(mipanalyzerobject_samples,
                                      mipanalyzerobject_gendistmat,
                                      type = c("mean", "median")
){

  if (is.list(mipanalyzerobject_gendistmat)) {
    gen_measure <- names(mipanalyzerobject_gendistmat)
    mipanalyzerobject_gendistmat <- mipanalyzerobject_gendistmat[[1]]
  } else {
    gen_measure <- deparse(substitute(mipanalyzerobject_gendistmat)) %>%
      strsplit("$", fixed = TRUE) %>%
      lapply(function(x) tail(x,1)) %>%
      unlist
  }

  smpl_admin_key1 <- mipanalyzerobject_samples %>%
    dplyr::select(c("sh312", "ADM1NAME")) %>%
    dplyr::mutate( item1 = as.numeric( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(smpl1 = sh312,
                  admin1_1 = ADM1NAME)


  smpl_admin_key2 <- mipanalyzerobject_samples %>%
    dplyr::select(c("sh312", "ADM1NAME")) %>%
    dplyr::mutate( item2 = as.numeric( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(smpl2 = sh312,
                  admin1_2 = ADM1NAME)


  gendist.tidy <-
    broom::tidy(as.dist( t( mipanalyzerobject_gendistmat) ))

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_admin_key1, by = "item1")

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_admin_key2, by = "item2")


  ret <- gendist.tidy %>%
    dplyr::group_by(admin1_1, admin1_2) %>%
    dplyr::summarise(
      gen_distance = mean(distance, na.rm = T)
    ) %>%
    dplyr::mutate(z = qnorm(percent_rank(gen_distance))) %>%
    dplyr::mutate(variable = gen_measure)

  return(ret)

}

# summarise gentic distance by cluster
get_cluster_gendist_summary <- function(mipanalyzerobject_samples,
                                        mipanalyzerobject_gendistmat,
                                        select = c("sh312", "ADM1NAME", "hv001"),
                                        type = c("mean", "median")){

  if (is.list(mipanalyzerobject_gendistmat)) {
    gen_measure <- names(mipanalyzerobject_gendistmat)
    mipanalyzerobject_gendistmat <- mipanalyzerobject_gendistmat[[1]]
  } else {
    gen_measure <- deparse(substitute(mipanalyzerobject_gendistmat)) %>%
      strsplit("$", fixed = TRUE) %>%
      lapply(function(x) tail(x,1)) %>%
      unlist
  }

  smpl_clust_key1 <- mipanalyzerobject_samples %>%
    dplyr::select(select) %>%
    dplyr::mutate( item1 = as.numeric( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(smpl1 = sh312,
                  admin1_1 = ADM1NAME,
                  clust_1 = hv001)


  smpl_clust_key2 <- mipanalyzerobject_samples %>%
    dplyr::select(select) %>%
    dplyr::mutate( item2 = as.numeric( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(smpl2 = sh312,
                  admin1_2 = ADM1NAME,
                  clust_2 = hv001)


  gendist.tidy <-
    broom::tidy(as.dist( t( mipanalyzerobject_gendistmat) ))

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_clust_key1, by = "item1")

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_clust_key2, by = "item2")

  ret <- gendist.tidy %>%
    dplyr::group_by(clust_1, clust_2) %>%
    dplyr::summarise(
      gen_distance = mean(distance, na.rm = T)
    ) %>%
    dplyr::mutate(variable = gen_measure) %>%
    dplyr::mutate(z = qnorm(percent_rank(gen_distance)))

  return(ret)

}

# summarise gentic distance by sample including cluster as that's the level
# of information we can append to these samples
get_raw_gendist_summary <- function(mipanalyzerobject_samples,
                                    mipanalyzerobject_gendistmat,
                                    select = c("hhid","hvidx","hv001"),
                                    type = c("mean", "median")){

  if (is.list(mipanalyzerobject_gendistmat)) {
    gen_measure <- names(mipanalyzerobject_gendistmat)
    mipanalyzerobject_gendistmat <- mipanalyzerobject_gendistmat[[1]]
  } else {
    gen_measure <- deparse(substitute(mipanalyzerobject_gendistmat)) %>%
      strsplit("$", fixed = TRUE) %>%
      lapply(function(x) tail(x,1)) %>%
      unlist
  }

  smpl_key1 <- mipanalyzerobject_samples %>%
    dplyr::select(select) %>%
    dplyr::mutate( item1 = as.numeric( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(hhid_1 = hhid,
                  hvidx_1 = hvidx,
                  clust_1 = hv001)


  smpl_key2 <- mipanalyzerobject_samples %>%
    dplyr::select(select) %>%
    dplyr::mutate( item2 = as.numeric( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(hhid_2 = hhid,
                  hvidx_2 = hvidx,
                  clust_2 = hv001)


  gendist.tidy <-
    broom::tidy(as.dist( t( mipanalyzerobject_gendistmat) ))

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_key1, by = "item1")

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_key2, by = "item2")

  ret <- gendist.tidy %>%
    dplyr::group_by(clust_1, clust_2) %>%
    dplyr::summarise(
      gen_distance = mean(distance, na.rm = T)
    ) %>%
    dplyr::mutate(variable = gen_measure) %>%
    dplyr::mutate(z = qnorm(percent_rank(gen_distance)))

  gendist.tidy$variable <- gen_measure
  gendist.tidy$z <- qnorm(percent_rank(gendist.tidy$gen_distance))

  return(gendist.tidy)

}


add_distance_summaries <- function(mipanalyzer_object) {

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results <- list()
  }

  if (is.null(mipanalyzer_object$results)) {
    mipanalyzer_object$results$raw_gen_dist <- list()
  }

  distances <- names(mipanalyzer_object$gen_distances)

  for (i in distances) {
    mipanalyzer_object$results$raw_gen_dist[[i]] <-
      get_raw_summaries(
        mipanalyzerobject_samples = mipanalyzer_object$samples,
        mipanalyzerobject_gendistmat = mipanalyzer_object$gen_distances[i],
        type = "mean" # TODO fix matcharg
      )
  }

  return(mipanalyzer_object)

}
