#########################################################################
# Purpose: Create a "Map" of starting parameters for the Cluster Inbreeding
#          Coefficient Gradient Descent algorithm to iterate over
#
# Author: Nicholas F. Brazeau
#
#########################################################################\
library(tidyverse)
#....................................................................................................
# SNAKEMAP for Clusters
#....................................................................................................
#..............................................................
# read in data
#..............................................................
# read in GE as import
mtdt <- readRDS("data/derived_data/sample_metadata.rds")
clsts <- sort(unique(c(mtdt$hv001)))
fs <- seq(0.1, 0.9, by = 0.1)
ms <- c(1e-12, 1e-10, 1e-8, 1e-6, 1e-5)

#....................................................................................................
# Start Parameters for clusters
#....................................................................................................
# cluster start params
fparams <- as.data.frame(matrix(fs, nrow = length(fs), ncol = length(clsts)))
colnames(fparams) <- clsts
fparams <- lapply(1:length(ms), function(x){return(fparams)}) %>%
  do.call("rbind.data.frame", .)
fparams <- fparams[order(fparams$`1`), ]
# all f and m params
clust_params <- cbind(fparams, "m" = ms)

#..............................................................
# Add in Learning Rate
#..............................................................
f_learningrate <- c(1e-7, 1e-6, 1e-5, 1e-4)
m_learningrate <- c(1e-18, 1e-15, 1e-12, 1e-10)
learningrates <- expand.grid(f_learningrate, m_learningrate)
colnames(learningrates) <- c("f_learningrate", "m_learningrate")

# cluster start params
clust_lrandparams <- lapply(1:nrow(learningrates), function(x){return(clust_params)}) %>%
  do.call("rbind.data.frame", .)
clust_lrandparams <- clust_lrandparams[order(clust_lrandparams$`1`), ]
clust_lrandparams <- cbind(clust_lrandparams, learningrates)

#..............................................................
# add in inputs
#..............................................................
lrandparams.gc <- clust_lrandparams
lrandparams.road <- clust_lrandparams
lrandparams.gc$inputpath <- "data/derived_data/clst_inbreeding_dat/gcdist_gens.RDS"
lrandparams.road$inputpath <- "data/derived_data/clst_inbreeding_dat/roaddist_gens.RDS"
lrandparams.gc$full_matrix <- FALSE
lrandparams.road$full_matrix <- FALSE
lrandparams <- rbind.data.frame(lrandparams.gc, lrandparams.road)

#..............................................................
# split these and write them out
#..............................................................
lrandparams.list <- split(lrandparams, 1:nrow(lrandparams))
paramsnames <- paste0("paramset_", 1:length(lrandparams.list), "_clust")

write_out_params <- function(paramset, outname, outdir){
  saveRDS(object = paramset, file = paste0(outdir, outname, ".RDS"))
}

outdir <- "data/derived_data/clst_inbreeding_dat/paramset_clust/"
dir.create(outdir, recursive = T)
mapply(write_out_params, paramset = lrandparams.list, outname = paramsnames, outdir = outdir)

#..............................................................
# make master map file
#..............................................................
snakemap <- data.frame(parampath = paramsnames)
mastermap <- cbind.data.frame(lrandparams, snakemap)
saveRDS(mastermap, "data/derived_data/clst_inbreeding_dat/paramset_clust/mastermap.RDS")
colnames(snakemap) <- c("#parampath")
write.table(x = snakemap,
            file = "data/derived_data/clst_inbreeding_dat/paramset_clust/snake_map.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)


#....................................................................................................
# SNAKEMAP for Voroni Province Tesselations
#....................................................................................................
#..............................................................
# read in data
#..............................................................
# read voroni tesselation territories as import
vrdf <- readRDS("data/distance_data/voroni_base.RDS")
provs <- sort(unique(c(vrdf$IPUMSID)))
fs <- seq(0.1, 0.9, by = 0.1)
ms <- c(1e-12, 1e-10, 1e-8, 1e-6, 1e-5)


# prov voroni tesselation start params
fparams <- as.data.frame(matrix(fs, nrow = length(fs), ncol = length(provs)))
colnames(fparams) <- provs
fparams <- lapply(1:length(ms), function(x){return(fparams)}) %>%
  do.call("rbind.data.frame", .)
fparams <- fparams[order(fparams$`1`), ]
# all f and m params
prov_params <- cbind(fparams, "m" = ms)

#..............................................................
# Add in Learning Rate
#..............................................................
f_learningrate <- c(1e-7, 1e-6, 1e-5, 1e-4)
m_learningrate <- c(1e-18, 1e-15, 1e-12, 1e-10)
learningrates <- expand.grid(f_learningrate, m_learningrate)
colnames(learningrates) <- c("f_learningrate", "m_learningrate")

# prov voroni tesselation start params
prov_lrandparams <- lapply(1:nrow(learningrates), function(x){return(prov_params)}) %>%
  do.call("rbind.data.frame", .)
prov_lrandparams <- prov_lrandparams[order(prov_lrandparams$`1`), ]
prov_lrandparams <- cbind(prov_lrandparams, learningrates)

#..............................................................
# add in inputs
#..............................................................
lrandparams.migrate <- prov_lrandparams
lrandparams.migrate$inputpath <- "data/derived_data/clst_inbreeding_dat/migrate_gens.RDS"
lrandparams.migrate$full_matrix <- TRUE
lrandparams <- lrandparams.migrate
#..............................................................
# split these and write them out
#..............................................................
lrandparams.list <- split(lrandparams, 1:nrow(lrandparams))
paramsnames <- paste0("paramset_", 1:length(lrandparams.list), "_prov")

write_out_params <- function(paramset, outname, outdir){
  saveRDS(object = paramset, file = paste0(outdir, outname, ".RDS"))
}

outdir <- "data/derived_data/clst_inbreeding_dat/paramset_prov/"
dir.create(outdir, recursive = T)
mapply(write_out_params, paramset = lrandparams.list, outname = paramsnames, outdir = outdir)

#..............................................................
# make master map file
#..............................................................
snakemap <- data.frame(parampath = paramsnames)
mastermap <- cbind.data.frame(lrandparams, snakemap)
saveRDS(mastermap, "data/derived_data/clst_inbreeding_dat/paramset_prov/mastermap.RDS")
colnames(snakemap) <- c("#parampath")
write.table(x = snakemap,
            file = "data/derived_data/clst_inbreeding_dat/paramset_prov/snake_map.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
