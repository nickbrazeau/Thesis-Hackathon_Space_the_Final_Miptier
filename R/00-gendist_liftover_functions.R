 getadmin_gendist_summary <- function(mipanalyzerobject_samples,
                                      mipanalyzerobject_gendistmat,
                                      type = c("mean", "median")){


  smpl_admin_key1 <- mipanalyzerobject_samples %>%
    dplyr::select(c("sh312", "ADM1NAME")) %>%
    dplyr::mutate( item1 = as.character( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(smpl1 = sh312,
                  admin1_1 = ADM1NAME)


  smpl_admin_key2 <- mipanalyzerobject_samples %>%
    dplyr::select(c("sh312", "ADM1NAME")) %>%
    dplyr::mutate( item2 = as.character( seq(1, nrow(.)) )) %>% # key
    dplyr::rename(smpl2 = sh312,
                  admin1_2 = ADM1NAME)


  gendist.tidy <-
    broom::tidy(as.dist( t( mipanalyzerobject_gendistmat) ))

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_admin_key1, by = "item1")

  gendist.tidy <- gendist.tidy %>%
    dplyr::left_join(x=., y =smpl_admin_key2, by = "item2")


  ret <- gendist.tidy %>%
    dplyr::group_by(admin1_1, admin1_2) %>%
    dplyr::summarise(
      adim_gen_dist = mean(distance, na.rm = T)
    )

  return(ret)

}
