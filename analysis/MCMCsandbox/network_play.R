# misc
clean_units <- function(x){
  attr(x,"units") <- NULL
  class(x) <- setdiff(class(x),"units")
  x
}


#..............................................................
# Dependencies
#..............................................................
library(igraph)
library(tidygraph)

# equation
# P(IBD_A) = F_A (exp(-d/M))

#..............................................................
# Setup
#..............................................................
# asssume Ne << N_ind
Ne <- 3 # unique markers at a given locus
nc <- 3 # number of clusters
lambda_i <- 3 # mean number of individuals per cluster
M <- 5e5 # migration rate scalar
F_in <- matrix(runif(nc^2, 0, 1),
              nrow = nc, ncol = nc)


#..............................................................
# Simulate Location and Prob
#..............................................................
# simulate locations
long <- runif(nc, min = -5, max = 5)
lat <- runif(nc, min = -5, max = 5)
coords <- sf::st_as_sf(data.frame(long = long, lat = lat), coords = c("long", "lat"), crs = 4326)
# calculate differences in distances
gc <- sf::st_distance(coords, which = "Great Circle")
# probabaility of distance
dist.prob <- clean_units( dexp(gc/M) )

# prob of pairwise IBD by cluster
clst.IBD <- F_in * dist.prob

#..............................................................
# Simulate Data
#..............................................................
# draw indiviauls per cluster
clst.ind <- rpois(nc, lambda_i) + 1 # fix so this is true ZTP


dat <- data.frame(
  clst = clst.ind.assn <- rep(1:nc, clst.ind),
  sampleid = seq(1:length(clst.ind))
)














