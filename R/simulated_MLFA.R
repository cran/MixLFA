#' simulated_MLFA: Simulated data from the MLFA model
#'
#' This dataset contains a list of matrices, each with a specific purpose or
#' structure. The list includes four matrices: \eqn{Y}, \eqn{X}, \eqn{Z}, and \eqn{id}. simulated using the following parameters (described in the simulation study in (Ounajim et al., 2023)):
#'
#' @format A list with four elements:
#' \describe{
#'   \item{Y}{the observed outcomes matrix with 10 columns and 500 rows (100 subjects with 5 observations each).}
#'   \item{X}{the fixed effects design matrix with two columns (two explanatory variables for explaining the factor loadings variation) and 500 rows.}
#'   \item{Z}{the random effect design matrix similar to X.}
#'   \item{id}{a vector of length 500 containing subject identifiers}
#' }
#' @usage data(simulated_MLFA)
#' @references
#' Ounajim, A., Slaoui, Y., Louis, P. Y., Billot, M., Frasca, D., & Rigoard, P. (2023). Mixture of longitudinal factor analyzers and their application to the assessment of chronic pain.
#' Statistics in medicine, 42(18), 3259–3282. https://doi.org/10.1002/sim.9804
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
"simulated_MLFA"
