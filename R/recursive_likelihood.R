#' @param distdata
#' @param scalar numeric; scalar for distance
#' @details distdata column names must be: (1) smpl1; (2) smpl2; (3) relatedness; (4) K1; (5) K2; (6) distance

# TODO take this into cpp so it is just distance matrix method
# confirm likelihood power/legitimacy

get_distance_geno_likelihood <- function(name, distdata, scalar, K1sub = NULL, K2sub = NULL){
  #..............................................................
  # Assertions
  #..............................................................
  # check that column names are named correctly
  # check that it is not an expanded distance matrix

  #..............................................................
  # Setup
  #..............................................................
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

  # find clusters for combinations which is different than clsts if we have subsetted
  if(is.null(K1sub) & is.null(K2sub)){
    f_ij.combns <- as.data.frame( t(combn(clsts, 2)) )
  } else {
    f_ij.combns <- as.data.frame( t(combn(c(K1sub, K2sub), 2)) )
  }


  f_ij <- apply(f_ij.combns, 1, function(x){
    ret <- mean( distdata$relatedness[c(distdata$K1 %in% as.vector(x) &
                                          distdata$K2 %in% as.vector(x) &
                                          distdata$K1 != distdata$K2)
                                      ] )
    return(ret)
  })

  #..............................................................
  # Get Migration Probability Matrix
  #..............................................................
  # here we still need a full migration matrix
  if(is.null(K1sub) & is.null(K2sub)){
    m_ij.combns <- as.data.frame( t(combn(clsts, 2)) )
  } else {
    m_ij.combns <- as.data.frame( t(combn(clsts, 2)) ) %>%
      dplyr::filter(c( V1 %in% c(K1sub, K2sub) | V2 %in% c(K1sub, K2sub) ))
  }

  m_dij <- apply(m_ij.combns, 1, function(x){
    ret <- unique( distdata$distance[c(distdata$K1 == x[1] &
                                         distdata$K2 == x[2] | # transitivity
                                         distdata$K2 == x[1] &
                                         distdata$K1 == x[2]
    )
    ] )
    return(ret)
  })


  if(is.null(K1sub) & is.null(K2sub)){
    #..............................................................
    # If full data, use distance matrices
    #..............................................................
    f_dij.dist <- matrix(NA, length(clsts), length(clsts))
    f_dij.dist[lower.tri(f_dij.dist)] <- f_ij # R work downs rows, so this will fill in appropriately
    # fill in the rest of the distance matrix
    diag(f_dij.dist) <- f_ii
    f_dij.dist[upper.tri(f_dij.dist)] <-  t(f_dij.dist[lower.tri(f_dij.dist)])



    m_dij.dist <- matrix(NA, length(clsts), length(clsts))
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
        if (i != j) {

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

          # clean up pij_final for out
          pij_final <- as.data.frame(pij_final)
          colnames(pij_final) <- clsts
          rownames(pij_final) <- clsts

          #..............................................................
          # Calculate log-likelihood
          #..............................................................
          LL <- sum(log(
            as.vector(pij_final[lower.tri(pij_final, diag = F)],
                      pij_final[upper.tri(pij_final, diag = F)]) # need both sides of matrix triangle since I have said that u != v, which means 1,2 and 2,1 are not the same
          ))

        } # end i
      } # end j
    }

  } else {
    #..............................................................
    # If subset data, use adjacency lists
    #..............................................................
    # set up
    f_ii <- f_ii[!is.nan(f_ii)]
    m_dij <- as.vector(unlist(m_dij))
    m_dij <- dexp(m_dij/scalar)
    m_dij.list <- split(m_dij, factor(m_ij.combns[,1]))
    # need to account for transitivity
    m_dij.list[[2]] <- append(m_dij.list[[2]], m_dij.list[[1]][[1]])


    for (i in 1:length(c(K1sub, K2sub))) {
      cat("Iteration i: ", i, "\n")
      for (j in 1:length(c(K1sub, K2sub))) {
        cat("     Iteration j: ", j, "\n")
        # first term
        p_ij1 <- 1* 1* f_ii # m_iu * m_ju * f_ii

        # second term, we have already made migration adj list to not include ii,
        # so can just loop through here
        p_ij2 <- 0
        for (v in  1:(length(clsts)-1)) {
          p_ijv <- 1 * m_dij.list[[j]][v] * f_ij # m_iu * m_jv * f_ij
          p_ij2 <- p_ij2 + p_ijv
        } # end V iter

        pij_final <- p_ij1 + p_ij2


      } # end i for adj list
    } # end j for adj list
    #..............................................................
    # Calculate log-likelihood
    #..............................................................
    LL <- log(pij_final)
  } # end if else

  #..............................................................
  # out
  #..............................................................
  ret <- list(
    pij = pij_final,
    LL = LL
  )
  return(ret)
}
