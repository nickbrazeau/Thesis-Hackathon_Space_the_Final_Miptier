
#------------------------------------------------
#' @title X-Y plot of genetic distance data vs spatial distance
#'
#' @description Plots for a given grouping variable all genetic distance vs
#'   spatial distance relationships for given bins
#'
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param select Character for the grouping selected
#' @param genetic_distances
#' @param spatial_distances
#' @param spatial_bins
#' @param genetic_bins
#' @param facet_x
#' @param facet_y
#'
#' @export

gen_spatial_dist_plot <- function(x,
                                  select = NULL,
                                  genetic_measures = NULL,
                                  spatial_measures = NULL,
                                  spatial_bins = NULL,
                                  genetic_bins = NULL,
                                  facet_x = "genetic_distances_variable",
                                  facet_y = "spatial_distances_variable") {

  # Checks
  # ----------------------------------------------------------------------------
  MIPanalyzer:::assert_custom_class(x, "mipanalyzer_biallelic")
  if (!is.null(spatial_bins)) {
    MIPanalyzer:::assert_numeric(spatial_bins)
  }

  if (!is.null(genetic_bins)) {
    MIPanalyzer:::assert_numeric(genetic_bins)
  }

  # check that distances exist
  if (!all(c("genetic_distances", "spatial_distances") %in% names(x$results))) {
    stop("Missing distance lists in mipanalyzer object results list")
  }

  MIPanalyzer:::assert_non_null(facet_x)
  MIPanalyzer:::assert_non_null(facet_y)

  # check that the measures are where we expect
  if (!is.null(genetic_measures)) {
    if (!all(genetic_measures %in% names(x$results$genetic_distances))) {
      stop(
        "Requested geneticmeasures are not found
         within x$results$genetic_measures"
        )
    }
  }

  if (!is.null(spatial_measures)) {
    if (!all(spatial_measures %in% names(x$results$spatial_distances))) {
      stop(
        "Requested spatial measures are not found
         within x$results$spatial_measures"
      )
    }
  }

  # Extra funcs
  # ----------------------------------------------------------------------------

  # helper function for binding lists of data frames
  rbind_list_base <- function(x) {
    x2 <- do.call(rbind.data.frame,
                  c(x, stringsAsFactors = FALSE, make.row.names = FALSE))
    rownames(x2) <- seq_len(dim(x2)[1])
    x2
  }

  # midpoints helper function
  midpoints <- function(a) {
    a[-length(a)] + diff(a) / 2
  }

  # Needed variable names
  # ----------------------------------------------------------------------------

  # create labels for plot
  xlab <- "Spatial Distance"
  ylab <- "Genetic Distance"

  # Combine our data
  # ----------------------------------------------------------------------------

  gs <- 1:(grep("genetic", x$results$genetic_distances[[1]] %>% names)[1]-1)
  gs <- names(x$results$genetic_distances[[1]])[gs]
  full_res <- dplyr::left_join(x = rbind_list_base(x$results$genetic_distances),
                               y = rbind_list_base(x$results$spatial_distances),
                               by = gs)

  # subset if needed
  if (!is.null(genetic_measures)) {
    full_res <-
      full_res[full_res$genetic_distances_variable %in% genetic_measures, ]
  }


  if (!is.null(spatial_measures)) {
    full_res <-
      full_res[full_res$spatial_distances_variable %in% spatial_measures, ]
  }


  # Bin data accordingly
  # ----------------------------------------------------------------------------

  # create our spatial breaks
  if (!is.null(spatial_bins)) {

    # what'sthe max spatial
    sp_max <- max(full_res$spatial_distances_mean, na.rm = TRUE)
    sp_breaks <- c(0, 1, seq(0, sp_max, length.out = spatial_bins - 1)[-1])

    # group by spatial
    full_res <- full_res %>%
      mutate(bin = cut(spatial_distances_mean,
                       sp_breaks,
                       include.lowest = TRUE))

    xlab <- "Binned Spatial Distance"
  }

  # create our genetic breaks
  if (!is.null(genetic_bins)) {
    gen_max <- max(full_res$genetic_distances_mean, na.rm = TRUE)
    gen_breaks <- c(0, seq(0, gen_max, length.out = spatial_bins)[-1])

    # group by genetics
    full_res <- full_res %>%
      mutate(bin = cut(genetic_distances_mean,
                       gen_breaks,
                       include.lowest = TRUE))
    ylab <- "Binnned Genetic Distance"
  }

  # Summarise data by our breaks and plot
  # ----------------------------------------------------------------------------

  # if we are not binning then just plot points
  if (is.null(genetic_bins) && is.null(spatial_bins) && is.null(select)) {

    plot <- ggplot(full_res, aes(x = spatial_distances_mean,
                                 y = genetic_distances_mean,
                                 weight = genetic_distances_length,
                                 size = genetic_distances_length)) +
      geom_point()

    # otherwise group accordingly
  } else {

    # first create our bin that also includes our selected variables
    if (!is.null(select)) {
      if (!c("bin" %in% names(full_res))) {
        full_res$bin <-
          interaction(full_res[, paste(select, 1:2, sep = "_")], sep = "|")
      }
    }

    nms_swap <- grep("[[:punct:]][[:digit:]]", names(full_res))
    l <- length(nms_swap)/2
    nms_swap <- c(nms_swap[c((l+1) : (2*l),1:l)], ((l*2)+1):length(names(full_res)))

    summarised_res <- data.table::rbindlist(
      list(
      full_res,
      full_res %>%
        set_names(names(full_res)[nms_swap]) %>%
        mutate(bin = interaction(full_res[, paste(select, 2:1, sep = "_")], sep = "|"))
      ),
      use.names = TRUE) %>%
      group_by_at(vars(bin,
                       !!facet_x,
                       !!facet_y,
                       genetic_distances_variable,
                       spatial_distances_variable)) %>%
      summarise(gen_dist = mean(genetic_distances_mean),
                gen_se = sd(genetic_distances_mean) / sqrt(n()),
                spat_dist = mean(spatial_distances_mean),
                spat_se = sd(spatial_distances_mean) / sqrt(n()),
                length = sum(genetic_distances_length)) %>%
      ungroup() %>%
      mutate(gen_UL = gen_dist + 1.96 * gen_se) %>%
      mutate(gen_LL = gen_dist - 1.96 * gen_se) %>%
      mutate(spat_UL = spat_dist + 1.96 * spat_se) %>%
      mutate(spat_LL = spat_dist - 1.96 * spat_se)

    plot <- ggplot(summarised_res, aes(x = spat_dist, y = gen_dist,size = length)) +
      geom_pointrange(aes(ymin = gen_LL, ymax = gen_UL)) +
      geom_errorbarh(aes(xmax = spat_UL, xmin = spat_LL, height = 0)) +
      geom_point(aes(size = length))

  }

  # add exponenetial fit and facet
  facet_x <- grep(facet_x, names(full_res), value=TRUE)
  facet_y <- grep(facet_y, names(full_res), value=TRUE)

  plot <- plot +
    geom_smooth(method="glm", aes(weight = length),
                formula = y+.Machine$double.xmin ~ x,
                method.args=list(family=gaussian(link="log"))) +
    ylab(ylab) +
    xlab(xlab) +
    theme_bw() +
    scale_size(name = "Size", range = c(0.1,0.5))

  if (length(unique(plot$data[[facet_y]])) == 1) {
    plot <- plot + facet_wrap(as.formula(paste("~", facet_x)), scales = "free_y") +
      ggtitle(paste0(facet_y, " = ", unique(plot$data[[facet_y]])))
  } else if (length(unique(plot$data[[facet_x]])) == 1 ) {
    plot <- plot + facet_wrap(as.formula(paste("~", facet_y))) +
      ggtitle(paste0(facet_x, " = ", unique(plot$data[[facet_x]])))
  } else {
    plot <- plot + facet_grid(as.formula(paste(facet_y ,"~", facet_x)),
                              scales = "free_y")
  }

  print(plot)
  invisible(plot)

}
