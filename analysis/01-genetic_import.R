#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import the genetic data used in the "big barcode" manuscript
# to use for the genetic-spatial analysis
#----------------------------------------------------------------------------------------------------
#.................................
# Dependencies
#.................................
gdrive <- tcltk::tk_choose.dir()

devtools::install_github("nickbrazeau/rplasmodium")
library(rplasmodium)
devtools::install_github("mrc-ide/MIPanalyzer")
library(MIPanalyzer)
library(vcfR)



#.................................
# read in filtered big-barcode panel
#.................................
mipbigpanel <- readRDS(paste0(gdrive, "/data/derived_data/biallelic_distances.rds"))
#.................................
# read in drug res panel
#.................................
mipDRpanel <- readRDS(paste0(gdrive, "...."))

biallelic_sites <- !stringr::str_detect(mipDRpanel$loci$ALT, ",")  # this works because we only have SNPs for now
# long way to subset to a biallelic from multiallelic vcf but checked in raw vcf on command line, this works fine for now
mipDRpanel$loci <- mipDRpanel$loci[biallelic_sites, ]
mipDRpanel$coverage <- mipDRpanel$coverage[,biallelic_sites]
mipDRpanel$counts <- mipDRpanel$counts[2,,biallelic_sites] # this is the wsnraf which is typically stored in the biallelic mip analyzer object
class(mipDRpanel) <- "mipanalyzer_biallelic"



#.................................
# error handle the overlapping sites
#.................................

if(any(
  paste0(mipbigpanel$loci[,1], mipbigpanel$loci[,2]) %in% paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2])
)){
  warning(paste("Overlapping loci in the big mip panel and DR panel -- there are", sum( paste0(mipbigpanel$loci[,1], mipbigpanel$loci[,2]) %in% paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2]) ),
                "overlapping sites between the DR panel and Big Panel VCF"))
}

mipBB_sub_repeats_loci <- which( paste0(mipbigpanel_sub$loci[,1], mipbigpanel_sub$loci[,2]) %in% paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2]) )
mipDR_sub_repeats_loci <- which( paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2]) %in%  paste0(mipbigpanel_sub$loci[,1], mipbigpanel_sub$loci[,2]) )

# merge the big barcode panel reads into DR reads
tempcounts <- mipbigpanel_sub$counts[, mipBB_sub_repeats_loci]
tempcounts[is.na(tempcounts)] <- 0
mipDRpanel$counts[, mipDR_sub_repeats_loci] <- mipDRpanel$counts[, mipDR_sub_repeats_loci] + tempcounts

tempcov <- mipbigpanel_sub$coverage[, mipBB_sub_repeats_loci]
tempcov[is.na(tempcov)] <- 0
mipDRpanel$coverage[, mipDR_sub_repeats_loci] <- mipDRpanel$coverage[, mipDR_sub_repeats_loci] + tempcov

# drop these sites in big barcode now
mipbigpanel_sub$counts <- mipbigpanel_sub$counts[, -c(mipBB_sub_repeats_loci)]
mipbigpanel_sub$coverage <- mipbigpanel_sub$coverage[, -c(mipBB_sub_repeats_loci)]
mipbigpanel_sub$loci <- mipbigpanel_sub$loci[-c(mipBB_sub_repeats_loci), ]





