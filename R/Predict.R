#' Predict the latent factors using the fixed effect coefficients from MLFA Results
#'
#' This function calculates the predicted latent factor conditional expectations given the covariates \eqn{X} for new observations. The predictions are obtained using the estimated fixed effect coefficients of the mixed effect model in the MLFA model.
#'
#' @param res_MLFA a list containing the MLFA model parameters returned by the MLFA function.
#' @param C an integer giving the number of mixture components.
#' @param d an integer giving the number of latent factors in the factor analysis model.
#' @param X a matrix containing the design matrix for fixed effects covariates for the mixed effect model explaining the latent factors (an unit column can be included to estimate an intercept). Default is the matrix used to estimate the model.
#' @import ggplot2
#' @export
#' @details
#' The function extracts the relevant coefficients from the \eqn{\beta} matrix for the specified class \eqn{C} and factor \eqn{d}.
#' These coefficients are then multiplied by the predictor matrix \eqn{X} to compute the predicted values.
#' The function does not perform any plotting; it simply returns the estimated conditional expectation of the latent factor \eqn{d}.
#'
#' @return A matrix containing the predicted values.
#'
#' @examples
#' # Load the necessary datasets
#' data(simulated_MLFA)  # Load a simulated dataset based on the MLFA model
#' # Extract matrices from the list
#' # Extract matrix Y of outcomes of interest for the factor analysis model
#' Y <- simulated_MLFA$Y
#' # Extract matrix X of fixed effect covariates for describing the latent factors
#' X <- simulated_MLFA$X
#' # Extract matrix Z of random effect covariates for describing the latent factors
#' Z <- simulated_MLFA$Z
#' # Extract matrix id containing subject identifiers.
#' id <-simulated_MLFA$id
#' #' # Run the MLFA (Mixture of Longitudinal Factor Analyzers) function with:
#' # C: number of classes or clusters in our simulated data was set to 2.
#' # d: number of latent factors in our simulated data was set to 1.
#' # max_it: maximum number of iterations is set to 50 for a quick test.
#' # Estimation of the parameters of the MLFA model using the simulated data.
#' result_MLFA <- MLFA(C = 2, d = 2, X, Y, Z, id, max_it = 50, fixed_factor =  c(1,6))
#' # predict the latent factor scores from the resulting MLFA model.
#' Predict(result_MLFA, 1, 1)

Predict <- function(res_MLFA, C, d , X=res_MLFA$X) {

  p = dim(X)[[2]]
  coefficientsC1 <- res_MLFA$beta[[C]][(1 + (d-1) * p):((d-1) * p + p), 1]
  result <- X  %*%  coefficientsC1

  return(result)
}
