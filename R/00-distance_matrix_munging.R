

# given a mipanalyzer samples and distance matrix create a i,j data.frame
# of distances, with the groupings determined by select and using the summary
# functions in funcs, with args to funcs passed by ...

# N.B. Note the args will be passed to all funcs... this should be fixed/tidied
get_dist_summary <- function(mipanalyzerobject_samples,
                             mipanalyzerobject_distmat,
                             dist_label = "genetic_distances",
                             select = c("ID"),
                             funcs = c("mean", "sd"),
                             ...) {
  # Checks
  # ----------------------------------------------------------------------------
  ## catch for incompatible length matches
  if (any(
    ncol(mipanalyzerobject_distmat) <
    apply(mipanalyzerobject_samples[, select, drop = FALSE], 2, function(x) {
      length(unique(x))
    })
  )) {
    stop(
      "Distance Matrix has fewer columns than unique
       groups selected for in the samples dataframe"
    )
  }

  # if we are passing the distance matrix in as a list then grab the name off
  # it otherwise grab the name via deparsing.
  # i.e. the name being whether this is a genetic or spatial distance matric
  if (is.list(mipanalyzerobject_distmat)) {
    measure <- names(mipanalyzerobject_distmat)
    mipanalyzerobject_distmat <- mipanalyzerobject_distmat[[1]]
  } else {
    measure <- deparse(substitute(mipanalyzerobject_distmat)) %>%
      strsplit("$", fixed = TRUE) %>%
      lapply(function(x)
        tail(x, 1)) %>%
      unlist
  }

  # checks on objects
  MIPanalyzer:::assert_dataframe(mipanalyzerobject_samples)
  MIPanalyzer:::assert_2d(mipanalyzerobject_distmat)
  MIPanalyzer:::assert_string(dist_label)
  MIPanalyzer:::assert_string(select)

  # further specific checks
  if (!all(select %in% names(mipanalyzerobject_samples))) {
    stop("selected variables not found within mipanalyzer sammples object")
  }

  # maybe more checks should be written to give better error messages on not
  # found functions or funtiona arguments

  # Create key data frames and handle labelling arguments
  # ----------------------------------------------------------------------------

  # create name for what type of distance measures we are dealing with
  measure_varname <- paste0(dist_label, "_variable")

  # create keys for the samples to be used for matching
  key1 <- mipanalyzerobject_samples %>%
    dplyr::select(select) %>%
    dplyr::mutate(item1 = as.numeric(seq(1, nrow(.))))

  # create names for the variable keys rquested in key 1
  nms1 <- paste0(select, "_", 1)
  names(key1)[names(key1) %in% select] <- nms1

  # create keys for the samples to be used for matching
  key2 <- mipanalyzerobject_samples %>%
    dplyr::select(select) %>%
    dplyr::mutate(item2 = as.numeric(seq(1, nrow(.))))

  # create names for the variable keys rquested in key 1
  nms2 <- paste0(select, "_", 2)
  names(key2)[names(key2) %in% select] <- nms2

  # turn our distance matric into a long data.frame and join by our keys
  gendist.tidy <-
    broom::tidy(as.dist(t(mipanalyzerobject_distmat)))

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x = ., y = key1, by = "item1")

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x = ., y = key2, by = "item2")

  # create vector of the names for our summary distance measures
  funcs <- c(funcs, "length")
  nms <- c(nms1, nms2,
           as.character(sapply(funcs, function(x) {
             paste0(dist_label, "_", x)
           })))

  # extra fuinction for calcuating z-score of quantiles
  z <- function(x) {
    qnorm(percent_rank(x))
  }

  # group only over the most unique select variable
  gr <- apply(
    mipanalyzerobject_samples[, select], 2, function(x) {
      length(unique(x))
    })


  # group by the requested vars, and apply each summary function
  # before renaming using nms and adding the measure variable

  # we do this either by realising that gendist.tidy is already the same
  # size as the result of this, i.e. we have chosen a select variable that is
  # at the individual level and just use gendist.tidy or we summarise

  if (c(max(gr)^2 / 2 - max(gr)/2) == nrow(gendist.tidy)) {

    # this is not ideal as i've just hard written this in.
    # but if there was a quick way to realise that because
    # the grouping used is at the individual level that any summary
    # stats are likely to be the same so we can generate it quickly.
    ret <- gendist.tidy[,nms[nms %in% names(gendist.tidy)]]
    ret[grep("mean", nms, value=TRUE)] <- gendist.tidy$distance
    ret[grep("sd", nms, value=TRUE)] <- NA
    ret[grep("length", nms, value=TRUE)] <- 1
    ret[measure_varname] <- measure

  } else {


    # group and apply our functions but only using the most
    # unique grouping rather than all as it's slow using
    # extra groups
    g1 <- nms1[which.max(gr)]
    g2 <- nms2[which.max(gr)]

    ret <- gendist.tidy %>%
      group_by_at(vars(c(g1, g2))) %>%
      summarize_at(vars(distance), .funs = sapply(funcs, match.fun), ...) %>%
      dplyr::ungroup()
    # dplyr::mutate_at(vars(names(.)[-c(1:(2 * length(select)))]),
    #                  .funs = list(z = ~ z(.))) %>%

    # work out which extra vars need to be added from gendist.tidy
    add <- c(seq_len(length(nms)))[-vapply(names(ret), grep, numeric(1), nms)]

    for(i in add) {
      nc <- nchar(nms[i])
      if (as.numeric(substr(nms[i], start = nc, stop = nc)) == 1){
        ret[nms[i]] <-
          gendist.tidy[match(ret[[g1]], gendist.tidy[[g1]]), nms[i]]
      } else {
        ret[nms[i]] <-
          gendist.tidy[match(ret[[g2]], gendist.tidy[[g2]]), nms[i]]
      }
    }

    ret <- ret[,order(vapply(names(ret), grep, numeric(1), nms))] %>%
      setNames(nms) %>%
      dplyr::mutate(!!measure_varname := measure)

  }

  return(ret)

}


# Loops through all distance matrices stored within the
# x (mipanalyzer_object) element supplied
# by `distances`, and then uses get_dist_summary
add_dist_summaries <- function(x,
                               select = "hv001",
                               distances = "genetic_distances",
                               funcs = c("mean", "sd"),
                               ...) {

  # Checks
  # ----------------------------------------------------------------------------
  MIPanalyzer:::assert_custom_class(x, "mipanalyzer_biallelic")
  MIPanalyzer:::assert_string(select)

  # is the distances object supplied the actual list or the name
  # and check accordingly
  if (inherits(distances, "list")) {
    d <- deparse(substitute(distances))
    distances <- tail(strsplit(d, "$", fixed = TRUE)[[1]], 1)
  }

  MIPanalyzer:::assert_string(distances)

  if (!c(distances %in% names(x))) {
    stop("Requested distances object not found in mipanalyzer object")
  }


  # Cretaing needed results lists and more checks on args
  # ----------------------------------------------------------------------------

  # do we need to make the results list
  if (is.null(x$results)) {
    x$results <- list()
  }

  # create the name of the list within results
  list_nm <- distances
  #list_nm <- paste0(select, "_", distances)

  # create this if needed
  if (is.null(x$results[[list_nm]])) {
    x$results[[list_nm]] <- list()
  }

  # what distance matrics are we looping over
  dist_mats <- names(x[[distances]])

  # create the summary for the distance matrices
  for (i in dist_mats) {
    x$results[[list_nm]][[i]] <-
      get_dist_summary(
        mipanalyzerobject_samples = x$samples,
        mipanalyzerobject_distmat = x[[distances]][i],
        select = select,
        dist_label = distances,
        funcs = funcs,
        ...
      )
  }

  return(x)

}




## These  two functions need to be better/slicker
## What seems sensible is to cache the function call with the
## returned object so that when we call append_genetic_distances it checks
## to see if we have already calculated the distance and thus does not need to
## But I run out of time to make it fancy like this.


append_genetic_distances <- function(x,
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

  # create the genetic_distances list if needed
  if (is.null(x$genetic_distances)) {
    x$genetic_distances <- list()
  }

  # Calculate distances
  # ----------------------------------------------------------------------------

  # calculate and store the results
  if (nm != "MLE_IBD") {
    x$genetic_distances[[nm]] <-
      dist_func(x, ...)
  } else {
    x$genetic_distances[[nm]] <-
      dist_func(x, ...)$mle
  }

  # attribute with function call (would be good to have the actual call)
  attr(x$genetic_distances[[nm]], "function_call") <-
    distance_functions[nm]

  return(x)

}


add_genetic_distances <- function(x,
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
    x <- append_genetic_distances(x = x,
                                  dist_func = distance_functions[i],
                                  unlist(args[i]))
  }

  return(x)

}
