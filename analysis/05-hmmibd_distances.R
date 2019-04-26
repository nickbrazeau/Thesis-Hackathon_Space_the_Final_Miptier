# hmmibd
# libraries and imports
gdrive <- tcltk::tk_choose.dir()
library(tidyverse)

# grab data and example couple of plots
drcmips <- readRDS(paste0(gdrive, "/data/derived_data/cd2013_gen_space_epi_final.rds"))

# read in hmmibd results
hmmibd_dist <- readRDS(paste0(gdrive,"/data/derived_data/hmmibd.rds"))
mipbigpanel <- readRDS(paste0(gdrive, "/data/derived_data/all_bigbarcodemips.rds"))

# turn into a distance matric
hmmibd_dist_mat <- matrix(NA, nrow = nrow(mipbigpanel$samples),
                          ncol = nrow(mipbigpanel$samples))

hmmibd_dist_viterbi_mat <- matrix(NA, nrow = nrow(mipbigpanel$samples),
                          ncol = nrow(mipbigpanel$samples))

count <- 1
count2 <- 0
for (i in 1:(nrow(mipbigpanel$samples)-1)) {

  count2 <- count2 + length(hmmibd_dist_mat[i,-(1:i)])
  hmmibd_dist_mat[i,-(1:i)] <- hmmibd_dist$fract_sites_IBD[count:count2]
  hmmibd_dist_viterbi_mat[i,-(1:i)] <- hmmibd_dist$fract_vit_sites_IBD[count:count2]
  count <- count2 + 1

}

# subset to DRC
drc <- mipbigpanel$samples$Country=="DRC"
hmmibd_dist_mat_drc <- hmmibd_dist_mat[drc,drc]
hmmibd_dist_viterbi_mat_drc <- hmmibd_dist_viterbi_mat[drc,drc]
drcmips$genetic_distances$HMMIBD <- hmmibd_dist_mat_drc
drcmips$genetic_distances$HMMIBD_viterbi <- hmmibd_dist_viterbi_mat_drc

# append it so it has the same format as the others
drcmips$results[["genetic_distances"]][["HMMIBD"]] <-
  get_dist_summary(
    mipanalyzerobject_samples = drcmips$samples,
    mipanalyzerobject_distmat = drcmips[["genetic_distances"]]["HMMIBD"],
    select = c("ADM1NAME","hv001","Barcode"),
    dist_label = "genetic_distances",
    funcs = c("mean", "sd")
  )

drcmips$results[["genetic_distances"]][["HMMIBD_viterbi"]] <-
  get_dist_summary(
    mipanalyzerobject_samples = drcmips$samples,
    mipanalyzerobject_distmat = drcmips[["genetic_distances"]]["HMMIBD_viterbi"],
    select = c("ADM1NAME","hv001","Barcode"),
    dist_label = "genetic_distances",
    funcs = c("mean", "sd")
  )

# look at highly related
highly_related <- drcmips$results$genetic_distances$MLE_IBD$genetic_distances_mean > 0.2

plot(drcmips$results$genetic_distances$MLE_IBD$genetic_distances_mean[highly_related],
     drcmips$results$genetic_distances$HMMIBD$genetic_distances_mean[highly_related])

plot(drcmips$results$genetic_distances$MLE_IBD$genetic_distances_mean,
     drcmips$results$genetic_distances$HMMIBD_viterbi$genetic_distances_mean)


##
l <- length(drcmips$genetic_distances$MLE_IBD %>% as.numeric %>% na.omit)
df <- data.frame("Distance"=c(drcmips$genetic_distances$MLE_IBD %>% as.numeric %>% na.omit,
                              drcmips$genetic_distances$HMMIBD %>% as.numeric %>% na.omit,
                              drcmips$genetic_distances$HMMIBD_viterbi %>% as.numeric %>% na.omit),
                 "Measure"=c(rep("MLE IBD",l), rep("HMMIBD",l), rep("HMMIBD_viterbi",l)))

g1 <- ggplot(df, aes(x=Distance, fill = Measure)) + geom_density(alpha=0.8) + ylab("Density") + xlim(c(0,0.125)) + scale_y_sqrt() + theme_bw()
g2 <- ggplot(df, aes(x=Distance, fill = Measure)) + geom_density(alpha=0.8) + ylab("Density") + scale_y_sqrt() + theme_bw()
gc <- cowplot::plot_grid(g1,g2)
