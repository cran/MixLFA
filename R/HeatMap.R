#' Generate Heatmap of an MLFA factor loadings
#'
#' This function generates a heatmap for visualizing the factor loadings from the MLFA model results.
#' @param res_MLFA a list containing the MLFA model parameters returned by the MLFA function.
#' @param C an integer Class to display
#' @import ggplot2
#' @import pheatmap
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
#' # Generate a heatmap of the factor loadings of the first cluster
#'generate_heatmap(result_MLFA,C=1)

generate_heatmap <- function(res_MLFA, C) {

  Lam1 <- res_MLFA$Lam[[C]]
  d <- dim(Lam1)[2]
  df <- data.frame()

  row_names <- rownames(Lam1)
  if (is.null(row_names)) {
    row_names <- as.character(1:nrow(Lam1))
  } else {
    row_names <- make.unique(row_names)
  }

  for (i in 1:d) {
    temp_df <- data.frame(
      value = Lam1[, i],
      variable = paste0("Factor", i),
      row = factor(row_names, levels = row_names)
    )
    df <- rbind(df, temp_df)
  }

  plot <- ggplot(df, aes(x = variable, y = row, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", high = "blue",mid = "white",
                        midpoint = 0) +
    labs(x = "Factor", y = "Variable", fill = "Loading") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
    print(plot)
}



