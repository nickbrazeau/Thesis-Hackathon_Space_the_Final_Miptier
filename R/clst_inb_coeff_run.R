#!/usr/bin/env Rscript

#..............................................................
# Purpose of this script is to run the cluster inbreeding
# coefficient gradient descent from MIPanalyzer
#
# Note, this function is not generalizable -- purpose is
# project specific
#..............................................................
library(discent)
library(optparse)

option_list=list(
  make_option(c("-m", "--mastermap"),
              type = "character",
              default = NULL,
              help = paste("start parameters path"),
              metavar = "character"),
  make_option(c("-s", "--seed"),
              type = "numeric",
              default = NULL,
              help = paste("Seed to use"),
              metavar = "numeric"),
  make_option(c("-l", "--mlowerbound"),
              type = "numeric",
              default = NULL,
              help = paste("Lower bound for M"),
              metavar = "numeric"),
  make_option(c("-u", "--mupperbound"),
              type = "numeric",
              default = NULL,
              help = paste("Upper bound for M"),
              metavar = "numeric"),
  make_option(c("-O", "--output"),
              type = "character",
              default = NULL,
              help = "Output directory path.",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$mastermap)){
  print_help(opt_parser)
  stop("Missing input argument", call. = FALSE)
}

if(is.null(opt$seed)){
  print_help(opt_parser)
  stop("Missing seed argument", call. = FALSE)
}

if(is.null(opt$mlowerbound)){
  print_help(opt_parser)
  stop("Missing lower bound argument for m", call. = FALSE)
}

if(is.null(opt$mupperbound)){
  print_help(opt_parser)
  stop("Missing upper bound argument for m", call. = FALSE)
}

if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Missing output argument", call. = FALSE)
}


set.seed(opt$seed)
#..............................................................
# extract and run
#..............................................................
params <- readRDS(file = opt$mastermap)
params <- unlist(params)
snake_start_params <- as.numeric( params[!names(params) %in% c("f_learningrate", "m_learningrate", "inputpath", "full_matrix")] )
names(snake_start_params) <- names(params[!names(params) %in% c("f_learningrate", "m_learningrate", "inputpath", "full_matrix")] )
input <- params["inputpath"]
input <- readRDS(input)
f_learningrate <- as.numeric( params["f_learningrate"] )
m_learningrate <- as.numeric( params["m_learningrate"] )
fullmatrix <- as.logical( params["full_matrix"] )

ret <- discent::deme_inbreeding_spcoef(K_gendist_geodist = input,
                                       start_params = snake_start_params,
                                       m_lowerbound = 1e-25,
                                       m_upperbound = 5e-4,
                                       f_learningrate = f_learningrate,
                                       m_learningrate = m_learningrate,
                                       full_matrix = fullmatrix,
                                       steps = 1e4,
                                       report_progress = TRUE)

#..............................................................
# out
#..............................................................
saveRDS(object = ret,
        file = opt$output)
