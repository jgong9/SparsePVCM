#' Fit a Sparse Principal Varying Coefficient Model
#'
#'\code{SparsePVCM} is a function to return the optimal SparsePVCM by choosing the one with the
#'smallest BIC value.
#'
#' @export SparsePVCM
#' @param y_sample The vector of a response variable
#' @param X_sample The design  matrix of covariates
#' @param U_sample The vector of an index variable
#' @param num_core The number of cores for parallel computing. If it is not given, the total number of cores minus 1 will be used.
#' @param K The number of cubic B-spline basis functions. If not given, the total number minus one will be used.
#' @param num_tau_1 The number of nonzero tuning parameter candidates for nuclear norm. Default is 6
#' @param num_tau_2 The number of nonzero tuning parameter candidates for group lasso. Default is 6
#' @param max_iter A numeric value for the maximum number of iterations
#' @param tol Tolerance value used as a stopping criterion
#' @return A list that contains the results of SPVCM estimation
#' @importFrom face select.knots
#' @importFrom splines splineDesign
#' @importFrom glmnet cv.glmnet
#' @importFrom parallel detectCores makeCluster parRapply stopCluster



SparsePVCM <- function(y_sample, X_sample, U_sample, K = 10, num_core = NULL,
                 num_tun_1 = 10, num_tun_2 = 10,
                 max_iter = 100, tol = 1e-5){

  if(is.null(num_core)){
    num_cores <- detectCores() -1
  } else if(num_core > detectCores()){
    stop("The given number of cores is greater than the actual cores in your machine.")
  }

  my_cluster <- makeCluster(num_cores)


  if(sum(is.integer(c(K, num_tun_1, num_tun_2, max_iter))) > 0 ){
    stop("Some arguments are not positive integers.")
  }

  n <- length(y_sample)
  p <- dim(X_sample)[2]
  if( is.null( colnames(X_sample)) ){
    colnames(X_sample) <- paste("X", 1:p, sep = "")
  }

  ## Step 1. Construct M_mat
  t_vec <- seq( min(U_sample), max(U_sample),length.out = 100)
  knots_t_vec <- select.knots(t_vec, knots = K - 3 ,p=3)
  B_mat <- (splineDesign(knots_t_vec, t_vec, ord = 4))

  G <- (t(B_mat) %*% B_mat ) / length(t_vec)
  G_sqrt <- eigen(G)$vectors %*% diag( sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  G_sqrt_inverse <- eigen(G)$vectors %*% diag( 1/sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  B_U_mat <- t( (splineDesign(knots = knots_t_vec, x = U_sample, ord = 4)) %*% G_sqrt_inverse)

  M_mat <- matrix(NA, n, K * p)
  for(i in 1:n){
    M_mat[i,] <- kronecker( t(B_U_mat)[i,] , X_sample[i,])
  }


  ## Step 2. Initial value by Ridge regression
  # 2-a. Tuning parameters for Theta_1
  temp <- cv.glmnet(x=M_mat, y=y_sample, alpha=0, intercept = F)
  theta_initial_vec <- coef(temp, s=temp$lambda.min)[-1]
  Theta_initial <- matrix(theta_initial_vec, p, K)
  svd_initial <- svd(Theta_initial)

  tau1_vec <- c( 0,
                 exp(seq(from = log(min(svd_initial$d)), to = log(svd_initial$d[1]), length.out = num_tun_1)) )

  col_norms_initial <- apply( (Theta_initial), 1, FUN = function(x){
    sqrt( sum( x^{2}  ) )
  } )
  tau2_vec <- c(0,
                exp(seq( log( min(col_norms_initial) ), log(max(col_norms_initial)) ,
                         length = num_tun_2) )
  )


  # 2-b. Tuning parameters and weights for Theta_2
  tau1_vec_2 <- c( 0,
                   exp(seq(from = log(min(svd_initial$d)^3), to = log(svd_initial$d[1]^3), length.out = num_tun_1)) )
  tau2_vec_2 <- c(0,
                  exp(seq( log( min(col_norms_initial)^2 ), log(max(col_norms_initial)^2) ,
                           length = num_tun_2) )
  )
  w1_vec <- svd_initial$d^(-2)
  w2_vec <- col_norms_initial^(-1)


  ## Step 3. Prepare for DYS
  # 3-a. step size of DYS
  MTM_eigen <- eigen( 1/n * t(M_mat) %*%M_mat)
  beta_val <- MTM_eigen$values[1]

  # 3-b. combinations of tuning parameters
  tuning = list()
  para_mat <- matrix(NA, length(tau1_vec) * length(tau2_vec), 2)


  len_tuning = 0
  for(l in tau1_vec){
    for(k in tau2_vec){
      len_tuning = len_tuning + 1
      tuning[[len_tuning]] = c(l, k)
      para_mat[len_tuning,] <- c(l, k)

    }
  }



  fit_list <- list()
  bic_vec <- vector("numeric", len_tuning)
  BIC_fit_list <- list()

  ## Step 4. Implement DYS
  fit_list <- parRapply(my_cluster, para_mat, FUN = dys_SparsePVCM ,
                        y = y_sample,
                        M = M_mat,
                        Theta_initial = Theta_initial,
                        rho = 1/beta_val,
                        w1 = rep(1, min(p,K)),
                        w2 = rep(1,p),

                        max_iter = max_iter, tol = tol )


  for(k in 1:len_tuning ){
    # fit_list[[k]] <- dys_spvcm(
    #   y = y_sample,
    #   M = M_mat,
    #   Theta_initial = Theta_initial,
    #   rho = 1/beta_val,
    #   tau = tuning[[k]],
    #   w1 = rep(1, min(p,K)),
    #   w2 = rep(1,p),
    #
    #   max_iter = max_iter, tol = tol )

    BIC_fit_list[[k]] <- bic_SparsePVCM(y=y_sample, M=M_mat, Theta_fit=fit_list[[k]])
    bic_vec[k] <- BIC_fit_list[[k]]$bic

  }

  ## Step 5. Choose the best model with bic
  loc_best_bic <-  which( min(bic_vec) == bic_vec )

  if(length(loc_best_bic) > 1){
    loc_best_bic <- loc_best_bic[length(loc_best_bic)] # Most penalized one
  }

  ## Step 6. Refinement process
  Theta_initial_2 <- fit_list[[loc_best_bic]]$Theta_final
  para_mat_2 <- matrix(NA, length(tau1_vec_2) * length(tau2_vec_2), 2)

  tuning_2 = list()
  len_tuning_2 = 0
  for(l in tau1_vec_2){
    for(k in tau2_vec_2){
      len_tuning_2 = len_tuning_2 + 1
      tuning_2[[len_tuning_2]] = c(l, k)
      para_mat_2[len_tuning_2,] <- c(l, k)

    }
  }

  fit_list_2 <- list()
  BIC_fit_list_2 <- list()
  bic_vec_2 <- vector("numeric", len_tuning_2)

  fit_list_2 <- parRapply(my_cluster, para_mat_2, FUN = dys_SparsePVCM ,
                          y = y_sample,
                          M = M_mat,
                          Theta_initial = Theta_initial_2,
                          rho = 1/beta_val,
                          w1 = w1_vec,
                          w2 = w2_vec,

                          max_iter = max_iter, tol = tol)

  for(k in 1:len_tuning_2 ){
    # fit_list_2[[k]] <- dys_spvcm(
      # y = y_sample,
      # M = M_mat,
      # Theta_initial = Theta_initial_2,
      # rho = 1/beta_val,
      # tau = tuning_2[[k]],
      # w1 = w1_vec,
      # w2 = w2_vec,
      #
      # max_iter = max_iter, tol = tol )

    BIC_fit_list_2[[k]] <- bic_SparsePVCM(y=y_sample, M = M_mat, Theta_fit=fit_list_2[[k]])
    bic_vec_2[k] <- BIC_fit_list_2[[k]]$bic
  }

  ## Step 7. Choose the best model with bic
  loc_best_bic_2 <-  which( min(bic_vec_2) == bic_vec_2 )
  if(length(loc_best_bic_2) > 1){
    loc_best_bic_2 <- loc_best_bic_2[length(loc_best_bic_2)] # Most penalized one
  }
  ## Final output

  output_list <- fit_list_2[[loc_best_bic_2]]
  X_selected <- colnames(X_sample)[ BIC_fit_list_2[[loc_best_bic_2]]$index_nonzero_row ]

  y_fitted <-   diag(X_sample %*% output_list$Theta_final %*% B_U_mat)

  stopCluster(my_cluster)


  return(
    list(
      "Theta" = output_list$Theta_final, "rank" = output_list$rank, "X_selected" = X_selected,
      "y_fitted"= y_fitted, "B_U_mat" = B_U_mat,
      "nonzero_index" = BIC_fit_list_2[[loc_best_bic_2]]$index_nonzero_row,
      "U" = output_list$U, "D" = output_list$D, "V" = output_list$V,
      "knots_t_vec"=knots_t_vec, "G_sqrt_inverse" = G_sqrt_inverse
    )
  )
}
