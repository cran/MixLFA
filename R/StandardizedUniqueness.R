#' Generate Standardized Uniqueness from MLFA Results
#'
#' This function generates uniqueness plots (proportion of variance in the outcome variables in Y explained by the factor analysis model) based on the estimated error variance.
#'
#' @param res_MLFA a list containing the MLFA model parameters returned by the MLFA function.
#' @param C an integer giving the number of mixture components.
#' @import ggplot2
#' @import GGally
#' @import dplyr
#' @export
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
#' # Generate the uniqueness plots for the first cluster
#' stdUnique(result_MLFA, C=1)
#'
#'

stdUnique <- function(res_MLFA, C) {

  # Vérification des éléments attendus dans res_MLFA
  if (!all(c("id", "ri", "Y", "tau") %in% names(res_MLFA))) {
    stop("res_MLFA must contain 'id', 'ri', 'Y', and 'tau'")
  }

  # Conversion de l'identifiant en vecteur si nécessaire
  id <- res_MLFA$id
  if (!is.vector(id)) {
    id <- as.vector(id)
  }

  # Assignation des classes en fonction des random intercepts
  classe <- apply(res_MLFA$ri, 1, FUN = which.max)

  unique_numbers <- unique(id)

  # Création d'une correspondance entre les identifiants et les classes
  class_mapping <- data.frame(number = unique_numbers, classe = classe)
  id_with_classes <- data.frame(number = id) %>%
    left_join(class_mapping, by = "number")

  Y2 <- data.frame(res_MLFA$Y)

  J <- dim(Y2)[[2]]

  # Vérification si Y2$class est bien rempli
  if (nrow(id_with_classes) != nrow(Y2)) {
    stop("The number of rows in id_with_classes does not match the number of rows in Y2")
  }

  # Ajout des colonnes de classe et d'identifiant patient
  Y2$class <- as.character(id_with_classes$classe)
  Y2$id_patient <- as.character(id_with_classes$number)

  # Vérification si la classe C existe dans Y2
  if (!any(Y2$class == C)) {
    stop(paste("Class", C, "not found in the data."))
  }

  # Calcul de la variance pour chaque variable mesurée dans la classe spécifiée
  variances_per_class <- apply(Y2[Y2$class == C, 1:J], FUN = var, MARGIN = 2)

  if (length(variances_per_class) != J) {
    stop("The length of variances_per_class does not match the number of variables (J).")
  }

  # Standardisation de la singularité
  phi <- data.frame(matrix(c(res_MLFA$tau[[C]]), byrow = TRUE, ncol = J, nrow = 1))

  if (nrow(phi) != 1 || ncol(phi) != J) {
    stop("Dimension mismatch: phi should have 1 row and J columns.")
  }

  phi[1, ] <- phi[1, 1:J] / variances_per_class

  # Génération du graphique des coordonnées parallèles
  p <- ggparcoord(data = phi, scale = "globalminmax",
                  columns = 1:J,
                  showPoints = TRUE, alphaLines = 0,
                  title = "Uniqueness of the measured variables for each class"
  ) +
    scale_color_manual(values = c("#00008b")) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 10)
    ) + ylab("Standardized uniqueness")

  print(p)
}
