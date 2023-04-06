#' Fitted values and estimated function evaluation
#'
#' A function to evaluate the estimated functions at some particular values of index variable and
#' to calculate the fitted values of response if new predictor values are given
#'
#' @export predict_SparsePVCM
#' @param fit The output of SparsePVCM
#' @param U_new A numeric vector of index variable values to evaluate estimated functions
#' @param X_new A design matrix of new predictor values that corresponds to \code{U_new}
#' @return A list that contains fitted response, evaluated coefficient functions and principal functions at given U_new
#' @importFrom splines splineDesign



predict_SparsePVCM <- function(fit, U_new, X_new=NULL){

  p <- dim(fit$Theta)[1]

  B_U_mat <- t( (splineDesign(knots = fit$knots_t_vec, x = U_new, ord = 4)) %*% fit$G_sqrt_inverse)
  gamma_mat <- fit$D %*% t(fit$V) %*% B_U_mat
  beta_mat <- fit$Theta %*% B_U_mat

  if(is.null(X_new)){
    y_new <- NULL
  } else {
    if(length(U_new) != dim(X_new)[1]){
      stop("The length of U_new is different from the number of rows of X_new.")
    }
    y_new <- diag( X_new %*% beta_mat)
  }
  return(list("y_new"=y_new, "beta_mat"=beta_mat, "gamma_mat"=gamma_mat))
}
