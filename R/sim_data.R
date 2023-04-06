#' Generate data for simulation studies
#'
#' @export sim_data
#' @param snr Signal-to-noise ratio
#' @return A list that contains all data with 3 different scenarios
#' @importFrom MASS mvrnorm



sim_data <- function(snr = 5){

for(jj in 1:3){
  for(ii in 1:3){
    seed_num <- 1234 + 10 * jj + 1 * ii
    set.seed(seed_num)
    ## seed setting for each different scenario

    ### Step 1. Parameter setting
    N <- 500 # the total number of datasets for each scenario
    n <- c(100,200,300)[jj]
    p <- c(8,16,24)[ii] # number of nonzero coefficient functions
    p0 <- c(8,16,24)[ii] # number of zero coefficient functions
    r_0 <- 4
    D_diagonal <- c(1,0.8,0.6,0.5) * sqrt(p)
    sigma2_true <- sum(D_diagonal^2) / (snr * p)

    cov_X_true <- matrix(NA, p+p0, p+p0)
    for(j1 in 1:(p+p0) ){
      for(j2 in 1:(p+p0) ){
        cov_X_true[j1, j2] <- (0.5)^abs(j1 - j2)
      }
    }

    A_mat_temp <- matrix(0, p, r_0)
    for(i in 1:r_0){
      A_mat_temp[(1+ (p/r_0) * (i-1) ):(p/r_0 * i), i] <- 1/ sqrt(c(2,4,6,12)[ii])
    }
    A_mat <- rbind(A_mat_temp, matrix(0,p0,r_0))
    D_mat <- diag(D_diagonal)


    ### Step 2. Data generation
    for(j in 1:N){
      U_sample <- runif(n, min= 0, max = 1)
      X_sample <- mvrnorm(n, mu=rep(0,p+p0), Sigma=cov_X_true )
      error_sample <- sqrt(sigma2_true) * rnorm(n,0, 1)

      for( i in 1:n){
        y_sample[i] <- X_sample[i,] %*% A_mat %*% D_mat %*%
          c(F_basis(1, U_sample[i]),F_basis(2, U_sample[i]),F_basis(3, U_sample[i]),F_basis(4, U_sample[i])) + error_sample[i]
      }

      data_list[[j]] <- data.frame(y_sample, X_sample, U_sample)
    } # total N iterations to generate N datasets
  }
}

  return(data_list)
}


