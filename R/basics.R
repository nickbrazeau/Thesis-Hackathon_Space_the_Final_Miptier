#' @title Logistic Transformation
logit <- function(x, tol=1e-4){
  return( log(((x+tol)/(1-x+tol))) )
}

#' @title Expit Transformation
expit <- function(x, tol=1e-4){
  return( 1/(1+exp(-x + tol)) )
}

#' @title scale as vector
my.scale <- function(x, ...){
  ret <- scale(x, ...)
  if(ncol(ret) == 1){
    return(as.numeric(ret)) # coerce matrix of 1 col to numeric vector
  } else{
    stop("Matrix has more than one column")
  }
}

#' @title Pretty Summary Table
make_pretty_summ_table <- function(dat, sumvar, groupingvar = "") {
  dat %>%
    group_by_at(groupingvar) %>%
    dplyr::summarise_at(., .vars = sumvar, .funs = list(n = length,
                                                        min = min,
                                                        median = median,
                                                        mean = mean,
                                                        stdev = sd,
                                                        max = max)) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    DT::datatable(., extensions='Buttons',
                  options = list(
                    searching = T,
                    pageLength = 20,
                    dom = 'Bfrtip',
                    buttons = c('csv')))


}


#' @title Quick Jpeg
jpgsnapshot <- function(outpath, plot, type = "wide") {
  assert_in(type, c("long", "wide"))
  if (type == "long") {
    jpeg(outpath, width = 8, height = 11, units = "in", res = 500)
    plot(plot)
    graphics.off()
  } else if (type == "wide") {
    jpeg(filename = outpath, width = 11, height = 8, units = "in", res = 500)
    plot(plot)
    graphics.off()
  }
}



#
