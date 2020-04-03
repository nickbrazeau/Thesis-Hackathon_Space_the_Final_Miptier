#!/usr/bin/env Rscript

#..............................................................
# Purpose of this script is to run the cluster inbreeding
# coefficient gradient descent from MIPanalyzer
#
# Note, this funciton is not generalizable -- purpose is
# project specfici
#..............................................................
remotes::install_github("nickbrazeau/MIPanalyzer", ref = "238ccb6c47a6d9c8148122a46f0b3bfe3f2b24c3")
library(MIPanalyzer)
library(optparse)

option_list=list(
  make_option(c("-m", "--mastermap"),
              type = "character",
              default = NULL,
              help = paste("start parameters path"),
              metavar = "character"),
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

if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Missing output argument", call. = FALSE)
}


#..............................................................
# extract and run
#..............................................................
params <- readRDS(file = opt$mastermap)
params <- unlist(params)
snake_start_params <- as.numeric( params[!names(params) %in% c("learningrate", "inputpath")] )
names(snake_start_params) <- names(params[!names(params) %in% c("learningrate", "inputpath")] )
input <- params["inputpath"]
input <- readRDS(input)
snake_learn <- as.numeric( params["learningrate"] )


ret <- MIPanalyzer::cluster_inbreeding_coef(clst_gendist_geodist = input,
                                            start_params = snake_start_params,
                                            m_lowerbound = 0,
                                            m_upperbound = 1,
                                            learningrate = snake_learn,
                                            steps = 1e3,
                                            report_progress = TRUE)

#..............................................................
# out
#..............................................................
saveRDS(object = ret,
        file = opt$output)
