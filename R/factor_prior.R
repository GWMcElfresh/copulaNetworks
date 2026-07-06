#' Check whether a correlation matrix is consistent with a low-rank factor structure
#'
#' @param corMatrix Symmetric correlation matrix (d x d).
#' @param nFactors Number of factors to assess (default 1).
#' @return List with eigenvalues, proportion of variance explained, and a
#'   warning flag when the first `nFactors` explain less than 50%.
#' @export
CheckFactorStructure <- function(corMatrix, nFactors = 1L) {
  cor_matrix <- as.matrix(corMatrix)
  eigen_decomp <- eigen(cor_matrix, symmetric = TRUE)
  eigenvalues <- eigen_decomp$values
  prop_var <- eigenvalues / sum(eigenvalues)
  cum_var <- cumsum(prop_var)
  low_factor_warning <- cum_var[nFactors] < 0.5
  if (isTRUE(low_factor_warning)) {
    warning(
      "First ", nFactors, " factor(s) explain only ",
      round(100 * cum_var[nFactors], 1),
      "% of correlation variance - factor model may be inadequate.",
      call. = FALSE
    )
  }
  list(
    eigenvalues = eigenvalues,
    propVar = prop_var,
    cumVar = cum_var,
    nFactors = nFactors,
    lowFactorWarning = low_factor_warning
  )
}

#' Fit a factor copula prior on a large reference sample (Phase 1)
#'
#' Requires the optional **FactorCopula** package. Continuous margins only in v1;
#' supply covariate-adjusted residuals before calling.
#'
#' @param data Data frame containing node columns.
#' @param nodeCols Character vector of variable names.
#' @param nFactors Number of latent factors (1 or 2).
#' @param linkingCopula Character vector of linking copula families per variable
#'   (e.g. `"bvn"` for Gaussian).
#' @param nQuad Number of quadrature points for latent integration.
#' @return List with `factorFit`, `loadings`, `logLik`, `nFactors`,
#'   `linkingCopula`, `marginalSpec`, and `impliedCor`.
#' @export
FitFactorCopulaPrior <- function(data,
                                 nodeCols,
                                 nFactors = 1L,
                                 linkingCopula = "bvn",
                                 nQuad = 25L) {
  if (!requireNamespace("FactorCopula", quietly = TRUE)) {
    stop(
      "Package 'FactorCopula' is required for FitFactorCopulaPrior(). ",
      "Install with install.packages('FactorCopula').",
      call. = FALSE
    )
  }

  input_matrix <- as.matrix(data[, nodeCols, drop = FALSE])
  if (!is.numeric(input_matrix)) {
    stop("FitFactorCopulaPrior() requires numeric continuous columns only.", call. = FALSE)
  }
  n_obs <- nrow(input_matrix)
  n_vars <- ncol(input_matrix)
  if (n_obs < 3 * n_vars) {
    warning(
      "Prior sample size (n = ", n_obs, ") is small relative to dimension (d = ", n_vars,
      "); factor estimation may be unstable.",
      call. = FALSE
    )
  }

  marginal_spec <- FitMarginalSpec(input_matrix, nodeCols = nodeCols)
  uniform_matrix <- UniformMarginalTransform(input_matrix)
  continuous_data <- qnorm(uniform_matrix)

  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop(
      "Package 'statmod' is required for FitFactorCopulaPrior(). ",
      "Install with install.packages('statmod').",
      call. = FALSE
    )
  }
  gl <- statmod::gauss.quad.prob(nQuad)
  cop_f1 <- rep(linkingCopula, n_vars)

  factor_fit <- if (nFactors == 1L) {
    FactorCopula::mle1factor(
      continuous = continuous_data,
      ordinal = NULL,
      count = NULL,
      copF1 = cop_f1,
      gl = gl,
      hessian = FALSE,
      print.level = 0
    )
  } else if (nFactors == 2L) {
    cop_f2 <- rep(linkingCopula, n_vars)
    FactorCopula::mle2factor(
      continuous = continuous_data,
      ordinal = NULL,
      count = NULL,
      copF1 = cop_f1,
      copF2 = cop_f2,
      gl = gl,
      hessian = FALSE,
      print.level = 0
    )
  } else {
    stop("nFactors must be 1 or 2 in v1.", call. = FALSE)
  }

  loadings <- factor_fit$taus
  names(loadings) <- nodeCols
  implied_cor <- ImpliedFactorCorrelation(loadings, nFactors = nFactors)

  list(
    factorFit = factor_fit,
    loadings = loadings,
    logLik = factor_fit$loglik,
    nFactors = nFactors,
    linkingCopula = linkingCopula,
    marginalSpec = marginal_spec,
    impliedCor = implied_cor,
    nodeCols = nodeCols,
    n = n_obs
  )
}

#' Build an approximate correlation matrix from factor loadings
#'
#' Converts Kendall tau loadings to Pearson correlations via the Gaussian
#' copula identity, then assembles a one- or two-factor correlation structure.
#'
#' @param loadings Named numeric vector of Kendall tau loadings per variable.
#' @param nFactors Number of factors used in fitting (1 or 2).
#' @return Symmetric correlation matrix.
#' @export
ImpliedFactorCorrelation <- function(loadings, nFactors = 1L) {
  lambda <- sin(pi * loadings / 2)
  n_vars <- length(loadings)
  var_names <- names(loadings)
  if (nFactors == 1L) {
    cor_mat <- outer(lambda, lambda)
    diag(cor_mat) <- 1
  } else {
    # ponytail: two-factor implied cor uses equal split across factors
    lambda_mat <- cbind(lambda / sqrt(2), lambda / sqrt(2))
    cor_mat <- lambda_mat %*% t(lambda_mat)
    diag(cor_mat) <- 1
  }
  colnames(cor_mat) <- rownames(cor_mat) <- var_names
  cor_mat
}

#' Simulate data from a fitted factor copula prior
#'
#' @param priorFit Output of [FitFactorCopulaPrior()].
#' @param nSim Number of simulated observations.
#' @param nObs Alias for `nSim` (either may be used).
#' @return Numeric matrix (nSim x d) of simulated continuous scores on the
#'   normal-copula scale.
#' @export
SimulateFactorCopula <- function(priorFit, nSim = NULL, nObs = NULL) {
  n_sim <- nSim %||% nObs
  if (is.null(n_sim)) {
    stop("Specify nSim or nObs.", call. = FALSE)
  }
  loadings <- priorFit$loadings
  lambda <- sin(pi * loadings / 2)
  n_vars <- length(loadings)
  f <- stats::rnorm(n_sim)
  eps <- matrix(stats::rnorm(n_sim * n_vars), nrow = n_sim, ncol = n_vars)
  sim_mat <- matrix(lambda, nrow = n_sim, ncol = n_vars, byrow = TRUE) * f
  sim_mat <- sim_mat + sweep(eps, 2, sqrt(pmax(1 - lambda^2, 1e-8)), `*`)
  colnames(sim_mat) <- names(loadings)
  sim_mat
}

#' Null-coalescing helper
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
