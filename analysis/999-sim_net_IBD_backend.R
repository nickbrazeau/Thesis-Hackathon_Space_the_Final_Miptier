# IBD Fit
# TODO should be doing a better job of sampling under distribution
# may need Bob's shape sampler to account for nonsymmetric high IBD chunks
xhat <- ibD$malecotf
ret <- MASS::fitdistr(xhat, "exponential")


# https://igraph.org/r/doc/erdos.renyi.game.html
# Generate random graphs according to the Erdos-Renyi model

n <- length(unique(c(ibD$smpl1, ibD$smpl2)))
edges <- choose(n, 2)
simnet <- igraph::erdos.renyi.game(n,
                                   p.or.m = edges,
                                   type = c("gnm"),
                                   directed = FALSE,
                                   loops = FALSE) %>%
  tidygraph::as_tbl_graph() %>%
  activate(nodes) %>%
  dplyr::mutate(names = 1:n)
