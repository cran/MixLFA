#' Extract Fixed Coefficients from MLFA Results
#'
#' This function extracts the fixed effect coefficients \eqn{\beta} from the results obtained from a
#' Mixture of Longitudinal Factor Analyzers (MLFA) model for a specified class and factor.
#'
#' @param res_MLFA list containing the MLFA model parameters returned by the MLFA function.
#' @param C an integer giving the number of mixture components.
#' @param d an integer giving the factor index from which to extract the coefficients. This corresponds to the specific latent factor of interest.
#' @export
#' @details
#' The function first determines the number of predictor variables (\eqn{p}) by evaluating the number of columns in the
#' predictor matrix \eqn{X} that was used in the MLFA. It then extracts the relevant coefficients from the estimated fixed effects \eqn{\beta} vector
#' associated with the specified class \eqn{C} and factor \eqn{d}. The \eqn{\beta} vector is structured such that coefficients
#' for each factor are stored in contiguous blocks; this function selects the appropriate block corresponding to the
#' factor \eqn{d} within the class \eqn{C}.
#'
#' @return A vector containing the coefficients for the specified class and factor.
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
#' # Extract the fixed effect coefficients for the latent factor 1 in cluster 1
#' coef_vector <- Fixed_coef(result_MLFA, C=1, d=1)
#'
Fixed_coef <- function(res_MLFA, C, d ) {

  # Get the number of predictors (p) from the dimension of X
  p = dim(res_MLFA$X)[[2]]

  # Extract the relevant coefficients for the specified class (C) and factor (d)
  result <- res_MLFA$beta[[C]][(1 + (d-1) * p):((d-1) * p + p), 1]

  # Return the extracted coefficient vector
  return(result)
}
