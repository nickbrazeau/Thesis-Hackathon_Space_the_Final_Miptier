#---------------------------------------------------------------------------
# Purpose of this script is to do fit our pairwise traditional
# epi models as a backend process WITHIN provinces.
# Don't want to do this every time we knit
#----------------------------------------------------------------------------
library(tidyverse)
library(MuMIn)
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


adm1.IBD.btwn <- ibD.long %>%
  dplyr::group_by(adm1name.x) %>%
  dplyr::summarise(
    n = n(),
    meanIBD = mean(malecotf),
    seIBD = sd(malecotf)/sqrt(n),
    LLIBD = meanIBD - 1.96 * seIBD,
    ULIBD = meanIBD + 1.96 * seIBD
  ) %>%
  dplyr::rename(adm1name = adm1name.x)

adm1.IBD.wthn <- ibD.long %>%
  dplyr::filter(adm1name.x == adm1name.y) %>%
  dplyr::group_by(adm1name.x) %>%
  dplyr::summarise(
    wthnn = n(),
    wthnIBD = mean(malecotf),
    sewthnIBD = sd(malecotf)/sqrt(wthnn),
    LLwthnIBD = wthnIBD - 1.96 * sewthnIBD,
    ULwthnIBD= wthnIBD + 1.96 * sewthnIBD
  ) %>%
  dplyr::rename(adm1name = adm1name.x)

adm1.IBD <- dplyr::left_join(adm1.IBD.btwn, adm1.IBD.wthn,
                             by = "adm1name")

#..............................................................
# Import Covariates for Province
#..............................................................
provCovar.raw <- readRDS("data/derived_data/covar_rasterstack_provlocations_raw.RDS")

# transform covars
provCovar <- provCovar.raw %>%
  dplyr::mutate(
    prev = my.scale(logit(prev, tol = 0.1)),
    precip = my.scale(precip),
    temp = my.scale(temp),
    elev = my.scale(elev),
    crops = my.scale(logit(crops, tol = 0.1)),
    actuse = my.scale(logit(actuse, tol = 0.1)),
    netuse = my.scale(logit(netuse, tol = 0.1)),
    housing = my.scale(logit(housing, tol = 0.1))
  )

adm1.IBD.covar <- dplyr::left_join(adm1.IBD, provCovar, by = "adm1name")

#..............................................................
# Fit Model and Perform Model Selection
#..............................................................
# https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet
# https://rstudio-pubs-static.s3.amazonaws.com/63556_e35cc7e2dfb54a5bb551f3fa4b3ec4ae.html
# https://bbolker.github.io/stat4c03/notes/mixed_details.html

covars <- c("prev", "precip", "temp",
            "elev", "crops", "actuse",
            "netuse", "housing")

# just make simple functions
adm1_models <- tibble::tibble(
  class = c("base", "base", "covars", "covars"),
  outcome = c("meanIBD", "wthnIBD", "meanIBD", "wthnIBD"),
  covars = list("", "", covars, covars),
  weights = c("n", "wthnn", "n", "wthnn")
)

# modeling equation
fit_glm_adm1 <- function(class, data, outcome, weights, covars){
  if (class == "base") {
    eq <- as.formula(paste0(outcome, "~ 1"))
  } else {
    eq <- as.formula(paste0(outcome, "~", paste(covars, collapse = "+")))
  }
  ret <- glm(eq,
             data = adm1.IBD.covar,
             family = gaussian(link="log"),
             na.action = "na.fail")
  return(ret)

}


#..............................................................
# Run Models
#..............................................................
adm1_models$fitout <- furrr::future_pmap(adm1_models, fit_glm_adm1)
adm1_models$modsum <- purrr::map(adm1_models$fitout, summary)


#..............................................................
# Explore Combinations of Saturated Models Models
#..............................................................
# limit combinations to 18, which is the number of our non-spatial covars times 2
adm1_models$dredging <- purrr::pmap(adm1_models[,c("class", "fitout")], function(class, fitout){
  if (class == "base") {
    ret <- NA
  } else {
    ret <- MuMIn::dredge(fitout, m.lim = c(1,18))
  }
  return(ret)
})

#..............................................................
# Save out
#..............................................................
dir.create("results/trad_epi_fits/", recursive = T)
saveRDS(adm1_models, "results/trad_epi_fits/Provglm_within_modfits.RDS")





