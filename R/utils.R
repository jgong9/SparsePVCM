##########################################
#####       Utility functions        #####
##########################################


### Fourier Basis
F_basis <- function(k, t){
  if(k==0){
    eval <- 1
  } else
    if (k==1){
      eval <- sqrt(2) * cos(2 * k * pi * t)
    } else if (k==2){
      eval <- sqrt(2) * sin(2 * (k-1) * pi * t)
    } else if (k==3){
      eval <- sqrt(2) * cos(2 * (k-1) * pi * t)

    } else if (k==4){
      eval <- sqrt(2) * sin(2 * (k-2) * pi * t)
    }
  return(eval)
}


### Inverse function to recover the original vector x from sort(x) and order(x)
inverse_order <- function(x_sorted, order){
  if(is.vector(x_sorted)){

    output <- vector("numeric", length(x_sorted))
    for (i in 1:length(x_sorted)){
      output[order[i]]<- x_sorted[i]
    }
    return(output)

  } else if(is.matrix(x_sorted)){
    output <- matrix(NA, dim(x_sorted)[1], dim(x_sorted)[2])
    for (i in 1:(dim(x_sorted)[1]) ){
      output[order[i],]<- x_sorted[i, ]
    }
    return(output)
  } else {
    stop("x_sorted must be either a vector or matrix.")
  }
}


### Function for parallel computing for VCM with bandwidth selection
kth_LOO_fun_row_par <- function(x, data_set = data_set, p = p, label_list,
                                h_candidate){

  K_h <- function(t, h){
    a <- t/h
    b <- 1-(a^2)
    b[ b<=0] <- 0
    res <- 0.75 * b / h

    return(res)
  }


  ## kth fold as a test set
  k <- label_list[[ as.numeric(x[1]) ]]
  data_train <- data_set[!k, ]
  data_test <- data_set[k, ]

  h_temp <- as.numeric(x[3])

  X_sample_temp <- data_train[,2:(1+p)]
  y_sample_temp <- data_train[,1]
  U_sample_temp <- data_train[,2+p]

  # Take evaluation point u
  u_vec = as.numeric(unique(data_set[,2+p]))
  ijk <- as.integer(x[2])
  u <- u_vec[ijk]

  # Return 0 if there is no u in U_test of the kth fold
  unique_U_test <- unique(data_test[,2+p])
  if( !(u %in% unique_U_test) ){
    return(0)
  } else {

    # Follow the linear equation of Fan and Zhang 2008
    U_u_mat <- diag(U_sample_temp - u )
    Gamma_u_mat <- cbind(X_sample_temp, U_u_mat %*% X_sample_temp)
    W_u_vec <- K_h(U_sample_temp-u, h=h_temp)

    fit_lm <- stats::lm.wfit(x = as.matrix(Gamma_u_mat), y = y_sample_temp, w = W_u_vec, singular.ok = TRUE)

    # Take the first p values for coefficient functions evaluated at u
    locLinear_beta <- fit_lm$coefficients[1:p]

    # Find y_i such that U_i == u in the test set
    loc_u_test<- data_test[,2+p] == u

    # Calculate a sum of square of  y_i - y_hat_i
    output <- sum( (data_test[loc_u_test,1] - data_test[loc_u_test, 2:(1+p)] %*% (locLinear_beta))^2)

    return(output)

  }

} # End of function for VCM with CV


### Epanechnikov kernel
K_h <- function(t, h){
  a <- t/h
  b <- 1-(a^2)
  b[ b<=0] <- 0
  res <- 0.75 * b / h

  return(res)
}
