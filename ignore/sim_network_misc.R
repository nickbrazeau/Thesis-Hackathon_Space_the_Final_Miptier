library(tidyverse)
library(tidygraph)
library(ggraph)



simnet <- play_smallworld(n_dim  = 1,
                dim_size = 10,
                order = 2,
                p_rewire = 0.05)

simnetplotobj <- simnet %>%
  ggraph(layout = 'kk') +
  geom_edge_link(width  = 1, color = "#969696") +
  geom_node_point(size = 4) +
  scale_color_continuous(guide = 'legend') +
  theme_graph()

plot(simnetplotobj)


# no weights
simnetplotobjCentrality <- simnet %>%
  activate(nodes) %>%
  dplyr::mutate(id = 1:10) %>%
  mutate(Centrality = centrality_pagerank()) %>%
  ggraph(layout = "graphopt") +
  geom_edge_link(color = "#000000", alpha = 0.2) +
  geom_node_point(aes(size = Centrality, colour = Centrality)) +
  geom_node_text(aes(label = id), repel = TRUE, size = 5)+
  scale_color_viridis_c("Centrality") +
  theme_graph()




# weights
nfb <- simnet %>%
  activate(edges) %>%
  as_tibble()
nfb$weight <- ifelse(nfb$from %in% c(5) | nfb$to %in% c(5),
                     1, NA)
nfb$weight[is.na(nfb$weight)] <- runif(20, min = 0.01, max = 0.3)


simnetplotobjCentrality2 <- simnet %>%
  activate(edges) %>%
  left_join(., nfb, by = c("from", "to")) %>%
  activate(nodes) %>%
  dplyr::mutate(id = 1:10) %>%
  mutate(Centrality = centrality_pagerank(weights = weight)) %>%
  ggraph(layout = "graphopt") +
  geom_edge_link(aes(width  = weight), color = "#000000", alpha = 0.2) +
  scale_edge_width("Weight", range = c(0.2, 1)) +
  geom_node_point(aes(size = Centrality, colour = Centrality)) +
  geom_node_text(aes(label = id), repel = TRUE, size = 5)+
  scale_color_viridis_c("Centrality") +
  theme_graph()


svglite::svglite(file = "~/Desktop/example_network_centrality_weights.svg",
                 width = 11, height = 7)

cowplot::plot_grid(simnetplotobjCentrality, simnetplotobjCentrality2, align = "v")
graphics.off()



# extract centrality
centrality <- simnet %>%
  tidygraph::activate(nodes) %>%
  dplyr::mutate(centrality = tidygraph::centrality_pagerank(weights = weight)) %>%
  tibble::as_tibble()



# community detection algorithms
simnet <- play_islands(n_islands = 5,
                       size_islands = 10,
                       p_within = 0.8,
                       m_between = 5)

simnetplotobjCDA <- simnet %>%
  activate(nodes) %>%
  dplyr::mutate(id = 1:50) %>%
  ggraph(layout = "kk") +
  geom_edge_link(color = "#000000", alpha = 0.2) +
  geom_node_point() +
  theme_graph()
plot(simnetplotobjCDA)


simnetplotobjCDA2 <- simnet %>%
  mutate(Community = as.factor(tidygraph::group_louvain())) %>%
  ggraph(layout = 'kk') +
  geom_edge_link(aes(alpha = ..index..), show.legend = FALSE) +
  geom_node_point(aes(colour = Community), size = 3) +
  theme_graph()

plot(simnetplotobjCDA2)

jpeg(file = "~/Desktop/example_network_community.jpg",
     width = 11, height = 7, units = "in", res = 500)

cowplot::plot_grid(simnetplotobjCDA, simnetplotobjCDA2, align = "v")
graphics.off()



