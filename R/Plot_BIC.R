#' Plot BIC Values Across Iterations
#'
#' This function plots the BIC values across iterations for convergence evaluation.
#'
#' @param res_MLFA a list containing the MLFA model parameters returned by the MLFA function.
#' @import ggplot2
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
#' # plot the BIC from iteration 1 to iteration \eqn{max_it}.
#' plot_BIC(result_MLFA)


plot_BIC <- function(res_MLFA) {

  all_BIC_values <- unlist(res_MLFA$BIC)

  iterations <- seq_along(all_BIC_values)

  BIC_data <- data.frame(Iteration = iterations, BIC = all_BIC_values)

  ggplot(BIC_data, aes(x = Iteration, y = BIC)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 2) +
    ggtitle("BIC Values Across Iterations") +
    xlab("Iteration") +
    ylab("BIC Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}
