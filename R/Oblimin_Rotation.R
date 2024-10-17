#' Apply Oblimin Rotation to MLFA Factor Loadings
#'
#'This function applies an oblimin rotation (from the package \code{GPArotation}) <doi:10.1177/0013164404272507> to the factor loadings from the results.
#' of a Mixture of Longitudinal Factor Analyzers (MLFA) model. The oblimin default parameters are used.
#'
#' @param res_MLFA a list containing the MLFA model parameters returned by the MLFA function.
#' @return A list similar to `res_MLFA`, but with the factor loadings rotated using the oblimin method.
#' @import GPArotation
#' @export
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
#' # Apply the oblimin rotation to the factor loadings
#' rotated_result <- Oblimin_Rotation(result_MLFA)

Oblimin_Rotation <- function(res_MLFA) {

  res_MLFA_final <- res_MLFA
  taille <- length(res_MLFA$Lam)

  for(i in 1:taille) {
    rotation_result <- oblimin(res_MLFA$Lam[[i]])
    res_MLFA_final$Lam[[i]] <- rotation_result$loadings
  }

  return(res_MLFA_final)
}


