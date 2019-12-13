#..............................................................
# Purpose of this script is to make permutation test for
# meiotic siblings with prevalence differences
#..............................................................
library(tidyverse)

#..............................................................
# Read in metadata
#..............................................................
drcsmpls <- readRDS("data/distance_data/drcsmpls_foruse.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("id", "hv001")) %>%
  dplyr::rename(name = id)

mtdt <- readRDS("data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(name = id) %>%
  dplyr::select(c("name", "country", "hv001", "adm1name", "longnum", "latnum")) %>%
  dplyr::filter(name %in% drcsmpls$name)
#..............................................................
# Read in IBD data
#..............................................................
ibD <- readRDS("data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD_polarized_biallelic_processed.long.rds")

# log2 IBD measure
ibD <- ibD %>%
  dplyr::mutate(malecotf_gens = -log2(malecotf),
                malecotf_gens_inv = 1/malecotf_gens)

#..............................................................
# Idenity samples with IBD that are greater than half the
# genome (i.e. at or above meiotic siblings)
# These represent likely recent transmission events
#..............................................................
ibD.meiotic <- ibD %>%
  dplyr::filter(malecotf_gens <= 1)
ibdmeioticsmpls <- as.character(ibD.meiotic$smpl1)
ibdmeioticsmpls <- c(ibdmeioticsmpls, as.character(ibD.meiotic$smpl2))





#..............................................................
# Permutation Testing
#..............................................................
# Wrapper for Meiotic Siblings
meiotic_sib_wrapper <- function(name, nsmpls, IBDdistrib, covardistrib){
  #' @param nsmpls numeric; number of samples to simulate as pairs for IBD
  #' @param IBDdistrib numeric vector; distribution of IBD in population
  #' @param covardistrib numeric vector; distribution of covariate in population

  IBDpermutator <- function(nsmpls, IBDdistrib, covardistrib){

    distmat <- matrix(NA, nrow = nsmpls, ncol = nsmpls)
    smpl.pairs <- sum(lower.tri(distmat, diag = F))
    # ibd
    smpl.pair.IBD <- sample(IBDdistrib, size = smpl.pairs, replace = T)
    # covar
    smpl.pair.covar <- tibble::tibble(
      covar.x = sample(covardistrib, size = smpl.pairs, replace = T),
      covar.y = sample(covardistrib, size = smpl.pairs, replace = T)
    )

    # convert to dist object
    distmat[lower.tri(distmat, diag = F)] <- smpl.pair.IBD
    # long out
    smpl.pair.IBD <- as.dist(distmat) %>%
      broom::tidy(.) %>%
      dplyr::bind_cols(., smpl.pair.covar) %>%
      magrittr::set_colnames(c("smpl1", "smpl2", "simibd", "covar.x", "covar.y"))
    return(smpl.pair.IBD)
  } # end nested function

  #..............................................................
  # Tidy up for meiotic siblings
  #..............................................................

  # given that you are at or above meiotic, this is what your
  # covar distribution should look like
  smpl.pair.IBD <- IBDpermutator(nsmpls, IBDdistrib, covardistrib)
  if (max(smpl.pair.IBD$simibd) >= 0.5) {
    smpl.pair.IBD.meiotic <- smpl.pair.IBD %>%
      dplyr::filter(simibd >= 0.5) %>%
      dplyr::mutate(covardiff = abs(covar.x) - abs(covar.y)) # make magnitude same always

    return(smpl.pair.IBD.meiotic$covardiff)

  } else {
    return(NA)
  }

} # end function


#..............................................................
# Bring in covariates for NULLs
#..............................................................
clstcovar <- readRDS("data/derived_data/covar_rasterstack_samplinglocations_raw.RDS")
prevdistrib <- clstcovar$prev
citydistrib <- clstcovar$urban
# find distribution of samples that are from same cluster
# subtraction here will be turned into a binary downstream (i.e. 0 means same cluster)
unique_clst_vs_same_distrib <- mtdt$hv001

#..............................................................
# Run Simulations
#..............................................................
simdf <- tibble::tibble(
  name = c("prev", "urban", "clstwthn"),
  nsmpls = nrow(drcsmpls),
  IBDdistrib = list(ibD$malecotf),
  covardistrib = list(prevdistrib, citydistrib, unique_clst_vs_same_distrib)
)
iters <- 1e4
simdf <- lapply(1:iters, function(x) return(simdf)) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(name)




# for slurm on LL
dir.create("results/meiotic_null_dist", recursive = T)
setwd("results/meiotic_null_dist/")
ntry <- 1028 # max number of nodes
sjob <- rslurm::slurm_apply(f = meiotic_sib_wrapper,
                            params = simdf,
                            jobname = 'meiotic_null',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = 4000,
                                                 array = sprintf("0-%d%%%d",
                                                                 ntry,
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1:00:00"))

cat("*************************** \n Submitted Permutation Sims \n *************************** ")







