#########################################################################
# Purpose: Collect the results from the cluster inbreeding calculation
#
# Author: Nicholas F. Brazeau
#
# Date: April 08 2020
#########################################################################
library(tidyverse)
#TODO account for which space we are using --> from input path
#..............................................................
# bring in mtdt for names
#..............................................................
clsts <- readRDS("data/derived_data/sample_metadata.rds") %>%
  dplyr::select(c("hv001")) %>%
  dplyr::filter(!duplicated(.))
#..............................................................
# bring in param master list
#.............................................................
mastermap <- readRDS("data/derived_data/clst_inbreeding_dat/paramset/mastermap.RDS")
mastermap.sm <- mastermap %>%
  dplyr::select(c("1", "m", "learningrate", "inputpath", "parampath")) %>%
  dplyr::rename(param_set = parampath,
                fstart = 1,
                mstart = m) %>%
  dplyr::mutate(
    distance = sub("_gens.RDS", "", sub("data/derived_data/clst_inbreeding_dat/", "", inputpath))
  ) %>%
  dplyr::select(-c("inputpath"))

#..............................................................
# read in data from slurm scr
#..............................................................
filepaths <- list.files("data/derived_data/clst_inbreeding_dat/clust_inbd_results/",
                        pattern = ".RDS", full.names = T)

read_results <- function(path, clstnames){
  dat <- readRDS(path)
  param_set_name <- gsub(".RDS", "", basename(path))

  ret <- as.data.frame(matrix(c(param_set_name, dat$Final_Fis, dat$Final_m),
                               nrow = 1))
  colnames(ret) <- c("param_set", clstnames, "m")
  return(ret)
}

clst_ret.list <- lapply(filepaths, read_results, clstnames = sort(clsts$hv001))

clst_ret <- clst_ret.list %>%
  dplyr::bind_rows()

mastermap.results <- dplyr::left_join(mastermap.sm, clst_ret, by = "param_set") %>%
  dplyr::select(c("fstart", "mstart", "learningrate", "param_set", "distance", "m",
                  dplyr::everything()))

#..............................................................
# explore results
#..............................................................
write_csv(mastermap.results, path = "~/Desktop/inb_clst_results.csv")

#..............................................................
# out
#..............................................................
write_rds(x = clst_ret, path = "data/derived_data/clst_inbreeding_dat/clst_inbreeding_runs.RDS")

#..............................................................
# play
#..............................................................
names(clst_ret.list) <- sapply(clst_ret.list, function(x){as.character(x$param_set)})

outtest <- readRDS(filepaths[grepl("paramset_199.RDS", filepaths)])
plot(outtest$m_run)
outtest.firun <- do.call("rbind.data.frame",  outtest$fi_run)
plot(outtest.firun[,3])

bettest <- readRDS(filepaths[grepl("paramset_216.RDS", filepaths)])
plot(bettest$m_run)
bettest.firun <- do.call("rbind.data.frame",  bettest$fi_run)
plot(bettest.firun[,3])


plot(outtest.firun[, 134])
plot(bettest.firun[, 134])
max(outtest.firun[99e2:1e4, 134] - lag(outtest.firun[9899:1e3, 134])[2:102])
max(outtest$m_run[99e2:1e4] - lag(outtest$m_run[9899:1e4])[2:102])

max(bettest.firun[99e2:1e4, 134] - lag(bettest.firun[9899:1e4, 134])[2:102])
max(bettest$m_run[99e2:1e4] - lag(bettest$m_run[9899:1e4])[2:102])

plot(outtest$m_run)
plot(bettest$m_run)


# sanity
