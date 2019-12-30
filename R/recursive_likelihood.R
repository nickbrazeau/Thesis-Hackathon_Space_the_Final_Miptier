#' @param distdata
#' @param scalar numeric; scalar for distance
#' @details distdata column names must be: (1) smpl1; (2) smpl2; (3) relatedness; (4) K1; (5) K2; (6) distance
get_distance_geno_likelihood <- function(distdata, scalar){
  #..............................................................
  # Assertions
  #..............................................................
  # check that column names are named correctly
  # check that it is not an expanded distance matrix

  #..............................................................
  # Setup
  #..............................................................
  source("R/pairwise_helpers.R")
  # Find unique clusters
  clsts <- sort(unique(c(distdata$K1, distdata$K2)))

  #..............................................................
  # Within cluster mean IBD Matrix
  #..............................................................
  # get mean fii
  f_ii <- sapply(clsts, function(x){
    ret <- mean( distdata$relatedness[c(distdata$K1 %in% x & distdata$K2 %in% x)] )
    return(ret)
  })
  f_ij.combns <- t(combn(clsts, 2))
  f_ij <- apply(f_ij.combns, 1, function(x){
    ret <- mean( distdata$relatedness[c(distdata$K1 %in% as.vector(x) &
                                          distdata$K2 %in% as.vector(x) &
                                          distdata$K1 != distdata$K2)
                                      ] )
    return(ret)
  })

  f_dij.dist <- matrix(NA, length(f_ij), length(f_ij))
  f_dij.dist[lower.tri(f_dij.dist)] <- f_ij # R work downs rows, so this will fill in appropriately
  # fill in the rest of the distance matrix
  diag(f_dij.dist) <- f_ii
  f_dij.dist[upper.tri(f_dij.dist)] <-  t(f_dij.dist[lower.tri(f_dij.dist)])


  #..............................................................
  # Get Migration Probability Matrix
  #..............................................................
  m_dij <- apply(f_ij.combns, 1, function(x){
    ret <- unique( distdata$distance[c(distdata$K1 == x[1] &
                                         distdata$K2 == x[2])
                                     ] )
    return(ret)
  })

  m_dij.dist <- matrix(NA, length(m_dij), length(m_dij))
  m_dij.dist[lower.tri(m_dij.dist)] <- m_dij # R work downs rows, so this will fill in appropriately
  # fill in the rest of the distance matrix
  diag(m_dij.dist) <- 0
  m_dij.dist[upper.tri(m_dij.dist)] <-  t(m_dij.dist[lower.tri(m_dij.dist)])
  m_dij.dist <- dexp(m_dij.dist/scalar)

  #..............................................................
  # Recursion
  #..............................................................
  pij_final <- matrix(NA, nrow = length(clsts), ncol = length(clsts)) # storage

  for (i in 1:length(clsts)) {
    cat("Iteration i: ", i, "\n")
    for (j in 1:length(clsts)) {
      cat("     Iteration j: ", j, "\n")

      # first term
      p_ij1 <- m_dij.dist[i,i] * m_dij.dist[j,i] * f_dij.dist[i,i]

      # second term
      p_ij2 <- 0
      for (v in  1:length(clsts)) {
        if(v == i){
          next
        } else {
          p_ijv <- m_dij.dist[i,i] * m_dij.dist[j,v] * f_dij.dist[i,j]
          p_ij2 <- p_ij2 + p_ijv
        }
      } # end V iter

      pij_final[i,j] <- p_ij1 + p_ij2

      #..............................................................
      # Calculate log-likelihood
      #..............................................................
      LL <- sum(log(
        as.vector(pij_final[lower.tri(pij_final, diag = T)])
      ))

    } # end i
  } # end j

  #..............................................................
  # out
  #..............................................................
  ret <- list(
    pij = pij_final,
    LL = LL
  )
  return(ret)
}
