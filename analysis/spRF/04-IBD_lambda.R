# find the lambda that describes the decay
# of IBD for each sample
# Will use the exponential distribution because it is
# consistent with theory and has one paramater, lambda
# that can desribe it
#............................

#............................
# Extract DRC Samples
#............................
DRCsmpl <- mtdt %>%
  dplyr::filter(country == "DRC")

# remove pairwise comparisons and just get smpls all IBD measures
ibd.temp <- ibd
colnames(ibd.temp) <- c("smpl", "smpl", "nsites", "malecotf")
ibd.DRC <- rbind(ibd.temp[,c(1,4)], ibd.temp[,c(2,4)]) %>%
  dplyr::filter(smpl %in% DRCsmpl$id)
ibd.DRC.list <- split(ibd.DRC, factor(ibd.DRC$smpl))




#............................
# Use Mass to get rate param
#............................
DRCsmpl.rateibd <- purrr::map(ibd.DRC.list, function(x){
  xhat <- x$malecotf * 100
  ret <- MASS::fitdistr(xhat, "exponential")
  return(ret)
})


#............................
# Extract Bits
#............................
DRC_ibddecay.df <- data.frame(id = levels(factor(ibd.DRC$smpl)),
                              rate = unlist( purrr::map(DRCsmpl.rateibd, "estimate") ) )



DRC_ibddecay.df <- left_join(x=DRCsmpl, y = DRC_ibddecay.df, by = "id") %>%
  dplyr::select(c("id", "hv001", "rate"))

ggplot() +
  geom_sf(data=DRCprov) +
  geom_point(data = DRC_ibddecay.df, aes(x=long, y=lat, color = rate), alpha = 0.4) +
  scale_color_viridis_c()


DRC_ibddecay.df %>%
  dplyr::left_join(., counts) %>%
  dplyr::filter(hv001 %in% c(1:50)) %>%
  ggplot() +
  geom_boxplot(aes(x=factor(hv001), y=rate, color = factor(n)))




#....................................................................................
# Wrangle Data for Spatial Random Forest
#....................................................................................
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>%
  magrittr::set_colnames(tolower(colnames((.)))) %>%
  dplyr::rename(hv001 = dhsclust) %>%
  dplyr::select(c("hv001", "geometry"))
ge.drcsmpls <- left_join(DRC_ibddecay.df, ge, by = "hv001")
ge.drcsmpls <- sf::st_as_sf(ge.drcsmpls)
great.circ.mat <-  geosphere::distm(x = sf::as_Spatial(ge.drcsmpls),
                                    fun = geosphere::distHaversine)

great.circ.df <- as.data.frame(great.circ.mat)
colnames(great.circ.df) <- paste0("gc", seq(1, ncol(great.circ.df)))

fm <- as.formula(paste("ibddecayrate ~ ", paste(colnames(great.circ.df), collapse = "+")))
DRC.sprf <- ranger::ranger(formula = fm,
                           data = cbind.data.frame(ibddecayrate = DRC_ibddecay.df$rate, great.circ.df),
                           quantreg=TRUE, num.trees=150, mtry = 20)



















