#' Fit a Varying Coefficient Model with a local polynomial regression
#'
#'\code{VCM} is a function to fit VCM with a local polynomial regression through Epanechnikov kernel.
#'
#' @export VCM
#' @param y_sample The vector of a response variable
#' @param X_sample The design  matrix of covariates
#' @param U_sample The vector of an index variable
#' @param index_i The integer vector of subject index ranging from 1 to the total number of subjects. Each entry corresponds to that of \code{y_sample}
#' @param U_new A vector of any numeric values within the range of \code{U_sample} to evaluate coefficient functions at the given values as \code{beta_Unew}
#' @param num_core The number of cores for parallel computing. If it is not given, the total number of cores minus 1 will be used.
#' @param num_fold The number of folds to implement K-fold cross validation to choose the optimal bandwidth
#' @param num_h The number of candidates for the optimal bandwidth computed based on Fan and Gijbels 1995
#' @return A list that contains the results of VCM estimation
#' @importFrom stats lm.wfit approx
#' @importFrom parallel detectCores makeCluster parRapply stopCluster


VCM <- function(y_sample, X_sample, U_sample, index_i, U_new = NULL,
                num_core=NULL, num_fold = 5, num_h=50 ){

  if(is.null(num_core)){
    num_cores <- detectCores() -1
  } else if(num_core > detectCores()){
    stop("The given number of cores is greater than the actual cores in your machine.")
  }

  my_cluster <- makeCluster(num_cores)


  n_total <- length(y_sample)
  test_index <- sort(unique(index_i), decreasing = F)
  n_subject <- length(test_index)
  p <- dim(X_sample)[2]

  if(length(index_i) != n_total){
    stop("y_sample and index_i must have the same length.")
  }

  if(sum( (1:n_subject) == test_index) != n_subject ){
    stop("index_i must have integers from 1 to the total number of subject.")
  }

  if(is.null(num_core)){
    num_cores <- detectCores() -1
  } else if(num_core > detectCores()){
    stop("The given number of cores is greater than the actual cores in your machine.")
  }

  ### 1. Sort data with respect to the increasing order of subject index
  h_candidate <- seq( (max(U_sample) - min(U_sample))/(n), (max(U_sample) - min(U_sample))/2, length.out = num_h)

  order_new <- sample(1:n_subject, n_subject,replace=F)

  ## Assign labels (1 to num_fold) for K-fold CV to subjects
  if(n_subject %%num_fold != 0 ){
    number_per_fold <- n_subject %/%num_fold
    add_last_fold <- n_subject %%num_fold

    label_vec <- rep(0, n_subject)
    for( k in 1:num_fold){
      label_vec[order_new[((k-1)*(number_per_fold) +1 ):( k * number_per_fold ) ]] <- k
    }

    label_vec[  order_new[  (num_fold * number_per_fold + 1):(num_fold * number_per_fold + add_last_fold  ) ] ] <- num_fold
  } else {
    number_per_fold <- n_subject %/%num_fold

    label_vec <- rep(0, n_subject)
    for( k in 1:num_fold){
      label_vec[order_new[((k-1)*(number_per_fold) +1 ):( k * number_per_fold ) ]] <- k
    }
  }

  ## Sort all data in order of index_i and
  ## assign the labels (for CV) to the individual value of y_sample
  order_vec <- order(index_i, decreasing = F)

  y_sample_sort <- y_sample[order_vec]
  X_sample_sort <- X_sample[order_vec, ]
  U_sample_sort <- U_sample[order_vec]

  freq_subject <- table(index_i[order_vec])
  label_vec_update <- rep(label_vec, times = as.integer(freq_subject))

  ### 2. Find the optimal bandwidth h
  ## Prepare for parallel computing
  data_set <- cbind(y_sample_sort,X_sample_sort, U_sample_sort)

  label_list_input <- list()
  num_Ui_for_test_folds <- vector("integer", length = num_fold)
  # it will be used to calculate MSE from SSE
  unique_U <- unique(U_sample_sort)

  for(iii in 1:num_fold){
    temp <- label_vec_update == iii
    label_list_input[[iii]] <- temp
    # num_Ui_for_test_folds[iii] <- sum(temp)
  }

  ## We NEED to evaluate beta(U_i)s for all i=1,\dots,num_total
  par_mat <- expand.grid(1:num_fold, 1:length(unique_U), h_candidate)


  ll_st <- Sys.time()
  aaa <- parRapply(my_cluster, par_mat, FUN = kth_LOO_fun_row_par , data_set = data_set, p = p, label_list =label_list_input)
  ll_end <- Sys.time()
  ll_end - ll_st

  SSE <- colSums( matrix(aaa, num_fold * length(unique_U), num_h) )

  # MSE <- SSE / num_Ui_for_test_folds
  best_ll <- which( (SSE)  == min(SSE)  )
  h_opt <- h_candidate[best_ll]

  ### 3. Fit the local polynomial with the optimal h
  ## Note that we changed the order of the orignal data
  locLinear_beta_mat_sort <- matrix(NA, n_total,p)
  ## u_vec <- seq(0,1,length.out= 18)

  for(jjj in 1:length(unique_U)){
    u <- unique_U[jjj]
    U_u_mat <- diag(U_sample_sort - u )
    Gamma_u_mat <- cbind(X_sample_sort, U_u_mat %*% X_sample_sort)
    W_u_vec <- K_h(U_sample_sort-u, h=h_opt)
    fit_lm <- lm.wfit(x = as.matrix(Gamma_u_mat), y = y_sample_sort, w = W_u_vec, singular.ok = TRUE)
    loc_u <- U_sample_sort == u
    locLinear_beta_mat_sort[ loc_u,] <- matrix( rep(fit_lm$coefficients[1:p], sum(loc_u)) , sum(loc_u), p, byrow = T)

    ## locLinear_beta_mat[ seq(jjj, n, by = 18),] <- matrix( rep(fit_lm$coefficients[1:p], dim(yeast$x)[1]) , dim(yeast$x)[1], p, byrow = T)
  }

  ### 4. Recover the original order
  locLinear_beta_mat <- inverse_order(locLinear_beta_mat_sort , order_vec)

  ## fitted value corresponding to y_sample
  y_fitted <- rowSums( X_sample * locLinear_beta_mat)

  ## R-squared
  R_squared <- 1- sum((y_sample - rowSums( X_sample * locLinear_beta_mat) )^2) /   sum(( y_sample  - mean(y_sample)  )^2)

  ## New betas(U_new) and y_hat for U_new and X_new based on linear approximation
  if(!is.null(U_new)){
    LL_coef_mat_approx <- matrix(NA, 100,p)
    for(iii in 1:p){
      lin_approx <- approx(x= U_sample, y= locLinear_beta_mat[,iii], xout= U_new, rule = 2 , ties = mean)
      LL_coef_mat_approx[,iii]<- lin_approx$y
    }
  } else {
    LL_coef_mat_approx <- NULL
  }

  stopCluster(my_cluster)


  return(list("y_fitted" = y_fitted, "beta_U" = locLinear_beta_mat, "beta_Unew" = LL_coef_mat_approx
              ))


} # End of VCM function
