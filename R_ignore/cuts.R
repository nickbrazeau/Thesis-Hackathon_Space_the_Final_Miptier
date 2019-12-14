### Long-Range Transmission Events and Urbanicity
```{r, results='asis'}

coords.pts.longrange <- coords.pts %>%
  dplyr::filter(hv001 %in% c(unique(long_trans$hv001.x), unique(long_trans$hv001.y)))

ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  #prettybasemap_nodrc_dark +
  geom_curve(data = long_trans, alpha = 0.5, size = 1.1,
             aes(x = long_jitter.x, y = lat_jitter.x,
                 xend = long_jitter.y, yend = lat_jitter.y,
                 color = malecotf_gens)) +
  geom_point(data = coords.pts.longrange, aes(x = long_jitter, y = lat_jitter,
                                              fill = urbanmean),
             size = 2, shape = 21, stroke = 0) +
  scale_color_viridis_c("IBD Generations", option="plasma", direction = -1) +
  scale_fill_distiller("Mean Urbanicity", type = "div", palette = "RdYlBu")

```

```{r, results='asis'}

coords.pts.longrange <- coords.pts %>%
  dplyr::filter(hv001 %in% c(unique(long_trans$hv001.x), unique(long_trans$hv001.y)))
urban <- raster::raster("~/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/urbanicity_raster/urbanicity.grd")

ggplot() +
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252", size = 0.05) +
  ggspatial::layer_spatial(data = urban, aes(fill = stat(band1))) +
  scale_fill_distiller("Urbanicity Score", type = "div", palette = "RdYlBu") +
  #prettybasemap_nodrc_dark +
  geom_curve(data = long_trans, alpha = 0.5, size = 1.1,
             aes(x = long_jitter.x, y = lat_jitter.x,
                 xend = long_jitter.y, yend = lat_jitter.y,
                 color = malecotf_gens)) +
  geom_point(data = coords.pts.longrange, aes(x = long_jitter, y = lat_jitter),
             size = 2, shape = 21, stroke = 0.1, fill = NA) +
  scale_color_viridis_c("IBD Generations", option="plasma", direction = -1)

```
