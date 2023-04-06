#' Davis-Yin splitting algorithm for SparsePVCM
#'
#'\code{dys_SparsePVCM} is a function to find the optimal solution with given
#'two weight vectors and two tuning parameters.
#'
#' @export dys_SparsePVCM
#' @param y The vector of a response variable
#' @param M Design matrix for the vectorized Theta matrix
#' @param Theta_initial A matrix for an initial value
#' @param w1 A numeric vector for the weights of nuclear norm
#' @param w2 A numeric vector for the adaptive weights of group lasso
#' @param rho The step size of DYS algorithm
#' @param tau A vector of two tuning parameters for nuclear norm and group lasso
#' @param max_iter A numeric value for the maximum number of iterations
#' @param tol Tolerance value used as a stopping criterion
#' @return A list that contains the results of DYS algorithm



dys_SparsePVCM <- function(y,  M, Theta_initial, w1, w2,
                     rho, tau, max_iter, tol){

  n <- length(y)
  p <- dim(Theta_initial)[1]
  K <- dim(Theta_initial)[2]
  n_iter <- 1
  converge <- 0
  dev_vec <- rep(NA, max_iter)
  rank_trace <- c()

  X_m <- Theta_initial
  Y_mplus1 <- Theta_initial
  ### DYS algorithm
  while(n_iter <= max_iter){

    ## Save (m) update of X
    Temp <- X_m

    ## (m+1) update of Y
    for(j in 1:p){
      x_j <- X_m[j,]
      x_j_norm <- sqrt( t(x_j) %*% x_j )

      if(tau[2]==0 ){
        y_j_update <-  x_j
      } else {
        if( x_j_norm <= (tau[2] * w2[j] * rho) ){
          y_j_update <- rep(0, K )
        } else {
          y_j_update <- c(1 - ( tau[2] * w2[j] * rho ) / x_j_norm ) * x_j

        }
      }

      Y_mplus1[j,] <- y_j_update
    }

    ## Prepare input of proximal operator
    delta_H_y_m <-  1/n * crossprod( M, M %*% c(Y_mplus1) - y )
    W_mat <- 2 * Y_mplus1 - rho * matrix(delta_H_y_m, p, K) - X_m

    ## (m+1) update of Z
    svd_W_mat <- svd(W_mat)
    temp_d <- (svd_W_mat$d - rho * tau[1] * w1)
    temp_d[ temp_d <= 0 ] <- 0
    rank <- sum( temp_d > 0 )

    Z_mplus1 <- svd_W_mat$u %*% diag(temp_d) %*% t(svd_W_mat$v)
    if( rank == 0){
      U <-NULL
      V <- NULL
      D <- NULL
    } else {
      U <- svd_W_mat$u[,1:rank]
      V <- svd_W_mat$v[,1:rank]
      D <- temp_d[1:rank]
    }

    ## (m+1) update of X
    X_mplus1 <- X_m +  Z_mplus1 - Y_mplus1
    X_m <- X_mplus1


    ## Stopping criteria
    deviation <- norm(X_m - Temp, "F") / max(1, norm(Temp, "F"))
    dev_vec[n_iter] <- deviation
    if(deviation <= tol){
      converge<-1
      break
    }
    n_iter <- n_iter + 1
    rank_trace <- c(rank_trace, rank)
  }
  ### Loop ends

  ## Exact sparse and low-rank Theta
  if(rank != 0){
    svd_sparse <- svd(Y_mplus1)

    U_final <- svd_sparse$u[,1:rank]
    D_final <- diag( x=svd_sparse$d[1:rank], nrow = rank, ncol = rank )
    V_final <- svd_sparse$v[,1:rank]

    Theta_final <- U_final %*% D_final %*% t(V_final)
  } else {
    Theta_final <- matrix(0, p,K)
    U_final <- NULL
    D_final <- NULL
    V_final <- NULL
  }

  return(
    list(
      "Theta_Y" = Y_mplus1, "Theta_Z" = Z_mplus1, "Theta_X" = X_mplus1,
      "Theta_final" = Theta_final,
      "rank" = rank, "U" = U_final, "V"=V_final, "D"= D_final, "converge"=converge, "rank_trace"=rank_trace,
      "iteration" = n_iter, "dev_vec"= dev_vec
    )
  )
}
