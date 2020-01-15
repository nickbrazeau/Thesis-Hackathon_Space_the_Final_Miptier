
clstcovar
urban1 <- clstcovar %>%
  dplyr::select(c("hv001", "urban")) %>%
  dplyr::rename(hv001.x = hv001)

urban2 <- clstcovar %>%
  dplyr::select(c("hv001", "urban")) %>%
  dplyr::rename(hv001.y = hv001)


highprop.pairs <- ibD %>%
  dplyr::left_join(., y = urban1, by = "hv001.x") %>%
  dplyr::left_join(., y = urban2, by = "hv001.y") %>%
  dplyr::filter(hv001.x != hv001.y) %>%
  dplyr::mutate(
    clstpair1 = ifelse(hv001.x < hv001.y, hv001.x, hv001.y),
    clstpair2 = ifelse(hv001.x < hv001.y, hv001.y, hv001.x),
    clstpair = paste0(clstpair1, "/", clstpair2)
  ) %>%
  dplyr::group_by(clstpair) %>%
  dplyr::summarise(
    n = n(),
    highprop = mean(malecotf >= 0.5),
    urbandiff = unique( (urban.x - urban.y)^2 ) # only one value since clst urb never changes
  )

highprop.pairs %>%
  dplyr::filter(highprop > 0 ) %>%
  ggplot() +
  geom_point(aes(x = urbandiff, y = highprop, size = n))





clsts <- unique( mtdt$hv001 )
clsts <- as.data.frame( t(combn(clsts, 2)) ) %>%
  magrittr::set_colnames(c("hv001.x", "hv001.y"))



clsts$prophighlyrelated <- purrr::pmap(clsts, function(hv001.x, hv001.y){
  dat <- ibD %>%
    dplyr::filter(hv001.x == hv001.x & hv001.y == hv001.y |
                    hv001.x == hv001.y & hv001.y == hv001.x) # transitivity
  ret <- mean(dat$malecotf >= 0.5)
  return(ret)
})
