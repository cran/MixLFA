#' @useDynLib MixLFA, .registration = TRUE
#' @import Rcpp
NULL

#' Estimates a Mixture of Longitudinal Factor Analyzers (MLFA) model
#'
#' This function performs a mixture of longitudinal factor analyzers on multivariate longitudinal data as described by Ounajim et al (2023). The MLFA model is written a two equations representing a measurement factor model describing the link between the outcomes / indicators of interest \eqn{Y} and one or several latent factors \eqn{eta}, and a structural mixed effect model describing the link between the latent factors and a set of explanatory variables in the design matrices \eqn{X} and \eqn{Z} for fixed and random effects respectively:
#' \deqn{y_{itj}=\sum_{c=1}^{C}\mathbb{1}_{\{v_{i}=c\}}\left(\Lambda_{jc} \eta_{i.tc}+\epsilon_{itjc}\right),}
#' \deqn{\eta_{iktc} = X_{iktc}\beta_{kc} +  Z_{iktc}b_{ikc} + e_{itc},}
#' where i is the subject index, t is the measurement time index, j is the outcome index and k is latent factor index. The model parameters are estimated using the Expectation-Maximization (EM) algorithm.
#'
#' @param C an integer giving the number of mixture components.
#' @param d an integer giving the number of latent factors in the factor analysis model.
#' @param Y a matrix containing the observed outcomes / indicators of interest for the factor model.
#' @param X a matrix containing the design matrix for fixed effects covariates for the mixed effect model explaining the latent factors (an unit column can be included to estimate an intercept).
#' @param Z a matrix containing the design matrix for random effects covariates for the mixed effect model explaining the latent factors.
#' @param id a vector containing subject identifiers.
#' @param fixed_factor a vector of integers of length d containing the columns in Y with factor loadings fixed to 1.
#' @param scale an optional Boolean indicating whether the matrix Y needs to be scaled or not.
#' @param max_it an integer giving the maximum number of iterations for the expectation-maximization algorithm for parameter estimation. The algorithm might stop before \eqn{max_it} if the mean absolute difference between two successive iterations is smaller than \eqn{10^(-5)}
#' @param seed a seed value for the random initialization of parameters.
#' @references
#' Ounajim, A., Slaoui, Y., Louis, P. Y., Billot, M., Frasca, D., & Rigoard, P. (2023). Mixture of longitudinal factor analyzers and their application to the assessment of chronic pain.
#' Statistics in medicine, 42(18), 3259–3282. https://doi.org/10.1002/sim.9804
#' @return a list with the following components:
#' \describe{
#'   \item{Lam}{a list of length \eqn{C} containing the estimated factor loading matrice for each cluster.}
#'   \item{beta}{a list of length \eqn{C} containing the estimated fixed effects coefficients vector of length p*d where p is the number of fixed effect and d is the number of latent factors. The function fixed_coef can be used to extract the fixed effect coefficient for a given cluster and for a given latent factor.}
#'   \item{S_b}{a list of length \eqn{C} containing the estimated covariance matrices of the random effects.}
#'   \item{S_e}{a list of length \eqn{C} containing the estimated covariance matrices of the error term in the mixed effect model (both intra and inter-factor covariances).}
#'   \item{tau}{a list of length \eqn{C} containing the vector of estimated variances of the error terms \eqn{\epsilon_{itjc}} in the factor 'measurement' model.}
#'   \item{pro}{a vector containing the estimated proportion of each cluster.}
#'   \item{BIC}{Bayesian Information Criterion for model selection for either the number of clusters or the number of latent factors.}
#'   \item{AIC}{Akaike Information Criterion for model selection.}
#'   \item{ICL}{Integrated Completed Likelihood for model selection.}
#'   \item{ri}{a matrix containig the the probability of class membership for each subject.}
#'   \item{VerifNan}{a Boolean indicating whether the model generated Nan values or not.}
#' }
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
#' # Print the resulting factor loading matrices and fixed effect coefficients
#' print(result_MLFA$Lam)
#' print(result_MLFA$beta)

MLFA <- function(C, d, X, Y, Z, id , max_it,fixed_factor=1:d, seed=1 , scale = TRUE ) {

  if( scale == TRUE)
  {
    Y <- scale(Y)
    warning("Y has been scaled")
  }

  fixed_factor <- fixed_factor - 1

  if (length(unique(fixed_factor)) != length(fixed_factor)) {
    stop("Error: The 'fixed_factor' parameter contains duplicate values. Please ensure all values are unique.")
  }

   if (length(fixed_factor) != d) {
    stop(paste("Error:fixed_factor must be of length", d, "but is of length", length(fixed_factor)))
  }

  if (any(fixed_factor < 0)) {
    stop("Error: All elements of fixed_factor must be positive or equal to 0.")
  }

  if (!is.matrix(X)) {
    stop("X must be a matrix.")
  }
  if (!is.matrix(Y)) {
    stop("Y must be a matrix.")
  }
  if (!is.matrix(Z)) {
    stop("Z must be a matrix.")
  }

  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.")
  }
  if (nrow(X) != nrow(Z)) {
    stop("X and Z must have the same number of rows.")
  }

  X <- cbind(id = id, X)
  Y <- cbind(id = id, Y)
  Z <- cbind(id = id, Z)

  X <- X[order(X[, "id"]), ]
  Y <- Y[order(Y[, "id"]), ]
  Z <- Z[order(Z[, "id"]), ]

  id<- id[order(id)]

  X <- data.matrix(X[, -1])
  Y <- data.matrix(Y[, -1])
  Z <- data.matrix(Z[, -1])

  pro_in = list()
  Lam_in = list()
  beta_in = list()
  S_b_in = list()
  S_e_in = list()
  tau_in = list()

  J=dim(Y)[[2]]
  p=dim(X)[[2]]
  n_subj=length(unique(id))

  if (any(fixed_factor >= J)) {
    stop(paste("Error: All elements of fixed_vector must be less than"," ",J,"."))
  }

  for (c in 1:C) {
    set.seed(c*seed)

    Lam_in[[c]] = matrix(rbinom(J * d, 1, 0.5), ncol = d, nrow = J)
    beta_in[[c]] = rnorm(p * d, 0, 1)
    Ss_b = matrix(runif(p^2 * d^2), nrow = p * d, ncol = p * d)

    S_b_in[[c]] = Ss_b %*% t(Ss_b)

    Ss_e = matrix(runif(d^2), nrow = d, ncol = d)

    S_e_in[[c]] = Ss_e %*% t(Ss_e)

    tau_in[[c]] = runif(J) + 2

    pro_in[[c]]  = rep(1 / c, c) + c(0.1, -0.1)

  }

  if (!is.list(beta_in) || !all(sapply(beta_in, is.matrix))) {
    beta_in <- lapply(beta_in, as.matrix)
  }

  if (!is.list(tau_in) || !all(sapply(tau_in, is.matrix))) {
    tau_in <- lapply(tau_in, as.matrix)

  }

  result<- .Call('_MixLFA_MLFA', C, d, X, Y, Z, id,  Lam_in, beta_in, tau_in, S_e_in, S_b_in, max_it, pro_in[[C]], n_subj, J, p,fixed_factor)
  for (c in 1:C) {
    rownames(result$Lam[[c]]) <- colnames(Y)

  }

  result <- append(result, list(Y = Y))
  result <- append(result, list(X = X))
  result <- append(result, list(id = id))

  return(result)
}









