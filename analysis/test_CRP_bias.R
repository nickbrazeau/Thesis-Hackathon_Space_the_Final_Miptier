
# test_CRP_bias.R
#
# Author: Bob Verity
# Date: 2020-01-30
#
# Purpose:
# To develop an intuition as to if/how estimating relatedness from a small
# sample of individuals in a wider network of related individuals has any
# systematic effect in terms of bias.
#
# Result:
# Estimates are unbiased. Sample size only affects precision.
#
# ------------------------------------------------------------------

# load bobfunctions2 package. If not already installed this can be obtained from
# Github via the command devtools::install_github('bobverity/bobfunctions2')
library(bobfunctions2)

# function to draw n values from Chinese restaurant process with concentration
# parameter theta
sim_crp <- function(n, theta) {
  x <- rep(NA, n)
  x[1] <- 1
  for (i in 2:n) {
    prob_seen_before <- tabulate(x)/(i-1+theta)
    prob_new <- theta/(i-1+theta)
    x[i] <- sample(1:(max(x, na.rm = TRUE)+1), 1, prob = c(prob_seen_before, prob_new))
  }
  ret <- as.dist(mapply(function(y) y == x, x))
  return(ret)
}

# define simulation parameters
n <- 2e2      # number of indivials in population
n_samp <- 50  # number of individuals to sample from population
f <- 0.6      # pairwise provability of being highly related (dictates concentration parameter through theta = 1/f-1)
reps <- 1e2   # number of times to repeat simulation

# initialise vectors for storing results at population and sample level
mean_IBD_pop <- mean_IBD_samp <- rep(NA, reps)
for (i in 1:reps) {
  sim1 <- sim_crp(n, 1/f-1)
  inds <- sample(1:n, n_samp)
  sim1_sub <- as.dist(as.matrix(sim1)[inds, inds])

  mean_IBD_pop[i] <- mean(sim1)
  mean_IBD_samp[i] <- mean(sim1_sub)
}

# plot population vs. sample
plot(mean_IBD_pop, mean_IBD_samp, xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, lty = 2)
