source("R/assertions_v5.R")

#' @title Make a nice DT table from a dataframe
#' @import dplyr, DT
pretty_DT_tab <- function(df) {
  df %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    DT::datatable(., extensions='Buttons',
                  options = list(
                    searching = T,
                    pageLength = 20,
                    dom = 'Bfrtip',
                    buttons = c('csv')))
}

#' @title Summarize a variable
#' @details Using dplyr scheme, summarize a specific variable
#' @param df dataframe; input
#' @param x character; name of summary var
#' @param groupingvar character; name of grouping var
#' @import dplyr

summ_var <- function(df, x, groupingvar) {
  p1 <- df %>%
    dplyr::group_by_at(groupingvar) %>%
    summarise_at(x, .funs = c(n = length, min = min, mean = mean, median = median, max = max))
  p2 <- df %>%
    dplyr::group_by_at(groupingvar) %>%
    summarise_at(x, .funs = c(LQ95 = quantile), probs = c(0.025))
  p3 <- df %>%
    dplyr::group_by_at(groupingvar) %>%
    summarise_at(x, .funs = c(UQ95 = quantile), probs = c(0.975))
  # out
  dplyr::left_join(p1, p2, by = groupingvar) %>%
    dplyr::left_join(., p3, by = groupingvar)
}



#' @title Logistic Transformation
#' @details  Convert positive numeric to logistic
#' @param x numeric vector
logit <- function(x, tol=1e-4){
  assert_bounded(x, left =  0, right = 1)
  assert_bounded(tol, left =  0, right = 1)
  # out
  return( log(((x+tol)/(1-x+tol))) )
}


#' @title Expit Transformation
#' @details Reverse the logistic transformation
#' @param x numeric vector
expit <- function(x, tol=1e-4){
  assert_numeric(x)
  assert_bounded(tol, left =  0, right = 1)
  # out
  return( 1/(1+exp(-x + tol)) )
}

#' @title Scale a variable
#' @details Scale a single variable and return as a numeric vector instead of a matrix
#' @param x numeric vector

my.scale <- function(x, ...){
  assert_numeric(x)
  ret <- scale(x, ...)
  if(ncol(ret) == 1){
    return(as.numeric(ret)) # coerce matrix of 1 col to numeric vector
  } else{
    stop("Matrix has more than one column")
  }
}


#' @title Quick Jpeg
jpgsnapshot <- function(outpath, plot, type = "wide") {
  assert_in(type, c("long", "wide"))
  if (type == "long") {
    jpeg(outpath, width = 8, height = 11, units = "in", res = 200)
    plot(plot)
    graphics.off()
  } else if (type == "wide") {
    jpeg(filename = outpath, width = 11, height = 8, units = "in", res = 200)
    plot(plot)
    graphics.off()
  }
}




