#########################################################################
# Purpose: Collect the results from the cluster inbreeding explorer
#
# Author: Nicholas F. Brazeau
#
# Date: April 08 2020
#########################################################################
library(tidyverse)

#..............................................................
# bring in param master list
#.............................................................

#..............................................................
# read in data from slurm scr
#..............................................................
filepaths <- list.files("~/Documents/MountPoints/mountedScratchLL/Projects/Space_the_Final_Miptier/clust_inbd_results/",
                        pattern = ".RDS", full.names = T)
read_results_quick <- function(path, clstnames){
  dat <- readRDS(path)
  param_set_name <- gsub(".RDS", "", gsub("paramset_", "", basename(path)))

  ret <- as.data.frame(matrix(c(param_set_name, dat$Final_Fis, dat$Final_m),
                               nrow = 1))
  colnames(ret) <- c("param_set", clstnames, "m")
  return(ret)
}

clst_ret.list <- lapply(filepaths, read_results_quick, clstnames = sort(clsts$hv001))

clst_ret <- clst_ret.list %>%
  dplyr::bind_rows()

#..............................................................
# out
#..............................................................
write_rds(x = clst_ret, path = "data/derived_data/clst_inbreeding_dat/clst_inbreeding_runs.RDS")




# sanity
