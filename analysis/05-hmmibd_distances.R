# hmmibd

# read in hmmibd results
hmmibd_dist <- readRDS(paste0(gdrive,"/data/derived_data/hmmibd.rds"))
mipbigpanel <- readRDS(paste0(gdrive, "/data/derived_data/all_bigbarcodemips.rds"))

# turn into a distance matric
hmmibd_dist_mat <- matrix(NA, nrow = nrow(mipbigpanel$samples),
                          ncol = nrow(mipbigpanel$samples))

count <- 1
count2 <- 0
for (i in 1:(nrow(mipbigpanel$samples)-1)) {

  count2 <- count2 + length(hmmibd_dist_mat[i,-(1:i)])
  hmmibd_dist_mat[i,-(1:i)] <- hmmibd_dist$fract_sites_IBD[count:count2]
  count <- count2 + 1

}

# subset to DRC
drc <- mipbigpanel$samples$Country=="DRC"
hmmibd_dist_mat_drc <- hmmibd_dist_mat[drc,drc]
drcmips$genetic_distances$HMMIBD <- hmmibd_dist_mat_drc

# append it so it has the same format as the others
drcmips$results[["genetic_distances"]][["HMMIBD"]] <-
  get_dist_summary(
    mipanalyzerobject_samples = drcmips$samples,
    mipanalyzerobject_distmat = drcmips[["genetic_distances"]]["HMMIBD"],
    select = c("ADM1NAME","hv001","Barcode"),
    dist_label = "genetic_distances",
    funcs = c("mean", "sd")
  )

# look at highly related
highly_related <- drcmips$results$genetic_distances$MLE_IBD$genetic_distances_mean > 0.2

plot(drcmips$results$genetic_distances$MLE_IBD$genetic_distances_mean,
     drcmips$results$genetic_distances$HMMIBD$genetic_distances_mean)
