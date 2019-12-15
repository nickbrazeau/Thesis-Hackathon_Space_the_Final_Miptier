#---------------------------------------------------------------------------
# Purpose of this script is to do fit our pairwise traditional
# epi models as a backend process. Don't want to do this every time we knit
#----------------------------------------------------------------------------
library(tidyverse)
library(MuMIn)
library(stargazer)
library(raster)
source("~/Documents/GitHub/Space_the_Final_Miptier/R/basics.R")
source("~/Documents/GitHub/Space_the_Final_Miptier/R/themes.R")
source("~/Documents/GitHub/Space_the_Final_Miptier/R/pairwise_helpers.R")
set.seed(48)

#..............................................................
# Import metadata
#..............................................................

drcsmpls <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/drcsmpls_foruse.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("id", "hv001")) %>%
  dplyr::rename(name = id)

mtdt <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/sample_metadata.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(name = id) %>%
  dplyr::select(c("name", "country", "hv001", "adm1name", "longnum", "latnum")) %>%
  dplyr::filter(name %in% drcsmpls$name)

# drc prov for plot
DRCprov <- sf::st_as_sf(readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/map_bases/gadm/gadm36_COD_1_sp.rds"))

#....................................................................................
# Import Genetic Data
#....................................................................................
ibD <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/bigbarcode_genetic_data/mipanalyzer.DRCibD_polarized_biallelic_processed.long.rds")

# log2 IBD measure
ibD <- ibD %>%
  dplyr::mutate(malecotf_gens = -log2(malecotf),
                malecotf_gens_inv = 1/malecotf_gens)
ibD.long <- expand_distance_matrix(ibD)


#..............................................................
# Wrangle outcomes
#..............................................................
# liftover IBD
mtdt.prov <- mtdt %>%
  dplyr::select(c("name", "adm1name"))

colnames(mtdt.prov)[1] <- "smpl1"
ibD.long <- dplyr::left_join(ibD.long, mtdt.prov, by = "smpl1")

colnames(mtdt.prov)[1] <- "smpl2"
ibD.long <- dplyr::left_join(ibD.long, mtdt.prov, by = "smpl2")


# loop through IBD now
adm1name.pairs <- tibble::as_tibble(
expand.grid(list(DRCprov$adm1name, DRCprov$adm1name)) ) %>%
  magrittr::set_colnames(c("adm1name.x", "adm1name.y")) %>%
  dplyr::filter(adm1name.x != adm1name.y)
# function for mean IBD btwn provs
get_mean_IBD_prov <- function(adm1name.x, adm1name.y){
  x = adm1name.x
  y = adm1name.y
  ret <- ibD.long %>%  dplyr::filter(adm1name.x == as.character(x) &
                                       adm1name.y == as.character(y)) %>%
    dplyr::summarise(
      n = n(),
      meanIBD = mean(malecotf),
      seIBD = sd(malecotf)/sqrt(n),
      LLIBD = meanIBD - 1.96 * seIBD,
      ULIBD = meanIBD + 1.96 * seIBD
    )
  return(ret)
}

adm1name.pairs$meanIBD <- purrr::pmap(adm1name.pairs, get_mean_IBD_prov)
adm1name.pairs <- adm1name.pairs %>%
  tidyr::unnest(cols = meanIBD)


#..............................................................
# Get Covariates for Provinces
#..............................................................
covarrstr <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/covar_rasterstack_raw.RDS")

extract_agg_raster_polygon <- function(rstrlyr, plygn){
  vals <- raster::extract(x = rstrlyr, y = sf::as_Spatial(plygn),
                          fun = mean,
                          na.rm = T,
                          sp = F
  )
  return(as.vector(vals))

}

adm1 <- DRCprov %>%
  dplyr::select(c("adm1name", "geometry"))
adm1.list <- split(adm1, 1:nrow(adm1))

provCovar <- lapply(adm1.list, extract_agg_raster_polygon, rstrlyr = covarrstr)
provCovar <- do.call("rbind.data.frame", provCovar)
colnames(provCovar) <- names(covarrstr)

#......................................
# Wrangle -- add back in for modeling
#......................................
provCovar <- cbind.data.frame(adm1, provCovar)
colnames(provCovar)[1] <- "adm1name.x"
adm1name.pairs <- dplyr::left_join(adm1name.pairs, provCovar, by = "adm1name.x")

colnames(provCovar)[1] <- "adm1name.y"
adm1name.pairs <- dplyr::left_join(adm1name.pairs, provCovar, by = "adm1name.y")


adm1name.pairs.diff <- adm1name.pairs %>%
  dplyr::mutate(
    prevdiff = (prev.x - prev.y)^2,
    precipdiff = (precip.x - precip.y)^2,
    tempdiff = (temp.x - temp.y)^2,
    elevdiff = (elev.x - elev.y)^2,
    cropsdiff = (crops.x - crops.y)^2,
    actusediff = (actuse.x - actuse.y)^2,
    netusediff = (netuse.x - netuse.y)^2,
    housingdiff = (housing.x - housing.y)^2
  ) %>%
  dplyr::select(c("adm1name.x", "adm1name.y",
                  "n", "meanIBD", "seIBD", "LLIBD", "ULIBD",
                  dplyr::ends_with("diff"))) %>%
  dplyr::mutate(
    # rescale everything
    prevdiff = my.scale(prevdiff),
    precipdiff = my.scale(precipdiff),
    tempdiff = my.scale(tempdiff),
    elevdiff = my.scale(elevdiff),
    cropsdiff = my.scale(cropsdiff),
    actusediff = my.scale(actusediff),
    netusediff = my.scale(netusediff),
    housingdiff = my.scale(housingdiff)
  )


#..............................................................
# Fit Model and Perform Model Selection
#..............................................................
# https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet
# https://rstudio-pubs-static.s3.amazonaws.com/63556_e35cc7e2dfb54a5bb551f3fa4b3ec4ae.html
# https://bbolker.github.io/stat4c03/notes/mixed_details.html

covars <- c("prevdiff", "precipdiff", "tempdiff",
            "elevdiff", "cropsdiff", "actusediff",
            "netusediff", "housingdiff")


mods <- tibble(
  class = rep(c("base", "full"), 3),
  space = rep(c("gcdistance.scale", "roaddistance.scale", "riverdistance.scale"), 2),
  covars = rep(c("", list(covars)), 3),
  outcome = "meanIBD"
) %>%
  dplyr::arrange(class)




#..............................................................
# Setup Data and Inputs
#..............................................................
space_prov_matrix <- readRDS("~/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/distancematrix_byprovince.rds") %>%
  dplyr::select(c(dplyr::starts_with("adm1name"), dplyr::everything()))


adm1name.pairs.diff <- long_distance_matrix_join(adm1name.pairs.diff, space_prov_matrix,
                                             by = c("adm1name.x", "adm1name.y"))

adm1name.pairs.diff <- adm1name.pairs.diff %>%
  dplyr::mutate(
                # note need to scale distances for stabilization of glmm
                gcdistance.scale = my.scale(gcdistance),
                roaddistance.scale = my.scale(roaddistance),
                riverdistance.scale = my.scale(riverdist)
                )




fit_glmm_adm1 <- function(class, outcome, covars, space){
  if (class == "base") {
    eq <- as.formula(paste0(outcome, "~",
                            "(1+", space, "|", "adm1name.x)"))
  } else {
    eq <- as.formula(paste0(outcome, "~", paste(covars, collapse = "+"),
                            "+ (1+", space, "|", "adm1name.x)"))
  }
  ret <- lme4::glmer(eq,
                     data =  adm1name.pairs.diff,
                     family = gaussian(link="log"),
                     na.action = "na.fail",
                     weights = n)
  return(ret)

}

#..............................................................
# Run Models
#..............................................................
mods$fitout <- furrr::future_pmap(mods, fit_glmm_adm1)
mods$modsum <- purrr::map(mods$fitout, summary)


#..............................................................
# Explore Combinations of Saturated Models Models
#..............................................................
# limit combinations to 18, which is the number of our non-spatial covars times 2
mods$dredging <- purrr::map(mods$fitout, function(x){
  if (length(x@beta) == 1) {
    ret <- NA
  } else {
    ret <- MuMIn::dredge(x, m.lim = c(1,18))
  }
  return(ret)
})

#..............................................................
# Save out
#..............................................................
dir.create("results/trad_epi_fits/", recursive = T)
saveRDS(mods, "results/trad_epi_fits/Provglmm_modfits.RDS")


