#' Map ranks to uniform pseudo-observations
#'
#' @param x Numeric vector.
#' @param n_ref Reference sample size for the denominator (default `length(x)`).
#' @return Numeric vector in (0, 1).
#' @keywords internal
rank_to_uniform <- function(x, n_ref = length(x)) {
  rank(x, ties.method = "average") / (n_ref + 1)
}

#' Transform each column to uniform margins via ranks
#'
#' @param inputMatrix Numeric matrix (n x p).
#' @param nRef Reference sample size for rank denominator (default nrow of matrix).
#' @return Matrix of pseudo-observations in (0, 1) with column names preserved.
#' @export
UniformMarginalTransform <- function(inputMatrix, nRef = NULL) {
  input_matrix <- as.matrix(inputMatrix)
  if (is.null(nRef)) {
    nRef <- nrow(input_matrix)
  }
  uniform_matrix <- apply(input_matrix, 2, function(x) {
    rank_to_uniform(x, n_ref = nRef)
  })
  colnames(uniform_matrix) <- colnames(input_matrix)
  uniform_matrix
}

#' Fit marginal specification from a prior (large) sample
#'
#' Stores empirical rank knots per column so update data can be mapped through
#' the prior ECDF. Covariate-adjusted marginals are deferred — supply
#' pre-adjusted residuals in `priorMatrix`.
#'
#' @param priorMatrix Numeric matrix or data frame (N x d) from the prior cohort.
#' @param nodeCols Character vector of column names (default: all columns).
#' @return List with `method`, `nRef`, `knots`, and `nodeCols`.
#' @export
FitMarginalSpec <- function(priorMatrix, nodeCols = colnames(priorMatrix)) {
  prior_matrix <- as.matrix(priorMatrix[, nodeCols, drop = FALSE])
  knots <- lapply(seq_len(ncol(prior_matrix)), function(j) {
    sort(prior_matrix[, j])
  })
  names(knots) <- colnames(prior_matrix)
  list(
    method = "empiricalRank",
    nRef = nrow(prior_matrix),
    knots = knots,
    nodeCols = colnames(prior_matrix)
  )
}

#' Map update values through a prior marginal specification
#'
#' @param updateMatrix Numeric matrix or data frame (n x d) from the update cohort.
#' @param marginalSpec Output of [FitMarginalSpec()].
#' @param nodeCols Character vector of columns to transform.
#' @return Matrix of pseudo-observations in (0, 1).
#' @export
ApplyMarginalSpec <- function(updateMatrix,
                              marginalSpec,
                              nodeCols = marginalSpec$nodeCols) {
  update_matrix <- as.matrix(updateMatrix[, nodeCols, drop = FALSE])
  n_ref <- marginalSpec$nRef
  uniform_list <- lapply(nodeCols, function(col) {
    apply_prior_ecdf(update_matrix[, col], marginalSpec$knots[[col]], n_ref)
  })
  uniform_matrix <- do.call(cbind, uniform_list)
  colnames(uniform_matrix) <- nodeCols
  uniform_matrix
}

#' Apply prior empirical CDF to a vector
#' @keywords internal
apply_prior_ecdf <- function(x, knots, n_ref) {
  vapply(x, function(v) {
    if (is.na(v)) {
      return(NA_real_)
    }
  # ponytail: empirical ECDF ceiling; upgrade path = parametric marginal CDFs
    (sum(knots < v, na.rm = TRUE) + 0.5 * sum(knots == v, na.rm = TRUE) + 0.5) / (n_ref + 1)
  }, numeric(1))
}
