#' Compute BIC value for SPVCM
#'
#'\code{bic_spvcm} computes the BIC value for a given SPVCM for model selection.
#'
#' @export bic_SparsePVCM
#' @param y The vector of a response variable
#' @param M Design matrix for the vectorized Theta matrix
#' @param Theta_fit The output list of dys_pvcm function
#' @return A list that contains BIC value and an index vector to locate nonzero coefficient functions

bic_SparsePVCM <- function(y, M, Theta_fit){

  n <- length(y)
  p <- dim(Theta_fit$Theta_final)[1]
  K <- dim(Theta_fit$Theta_final)[2]
  rank_est <- Theta_fit$rank

  index <- apply(Theta_fit$Theta_Y==0, 1, sum)
  index_nonzero_row <- (index != K)
  p0 <- min(sum(index_nonzero_row), p)

  df <- K * rank_est  + p0 * rank_est - rank_est * rank_est

  if(rank_est == 0 | p0 == 0){
    SSE <- sum( (y - mean(y))^2 )
    bic_val <-  SSE / (n) # follow Jiang's BIC (2013) for the purpose of completeness
  } else {
    SSE <- sum( (y - M %*% c(Theta_fit$Theta_final)  )^2 )
    bic_val <- log(SSE / (n)) + df * log(n) / (n)
  }
  return(list(bic=bic_val, index_nonzero_row = index_nonzero_row))
}
