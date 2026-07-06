#' Pairwise Kendall tau deviation from implied prior correlations
#'
#' @param uniformMatrix Matrix of uniform pseudo-observations (n x d).
#' @param impliedCor Prior implied correlation matrix.
#' @return Mean absolute deviation of sample taus from implied taus.
#' @keywords internal
pairwise_tau_stat <- function(uniformMatrix, impliedCor) {
  vars <- colnames(uniformMatrix)
  if (is.null(vars)) {
    vars <- paste0("V", seq_len(ncol(uniformMatrix)))
    colnames(uniformMatrix) <- vars
  }
  implied_cor <- impliedCor[vars, vars, drop = FALSE]
  implied_tau <- (2 / pi) * asin(implied_cor)
  diag(implied_tau) <- NA

  pairs <- utils::combn(vars, 2, simplify = FALSE)
  if (length(pairs) == 0) {
    return(0)
  }

  obs_tau <- vapply(pairs, function(pr) {
    stats::cor(uniformMatrix[, pr[1]], uniformMatrix[, pr[2]], method = "kendall")
  }, numeric(1))

  impl_tau <- vapply(pairs, function(pr) {
    implied_tau[pr[1], pr[2]]
  }, numeric(1))

  mean(abs(obs_tau - impl_tau), na.rm = TRUE)
}

#' Test whether update data are consistent with the Phase 1 factor prior
#'
#' Simulation-based goodness-of-fit: compares a pairwise Kendall tau statistic
#' on the real update sample to a null distribution from prior replicates.
#'
#' @param priorFit Output of [FitFactorCopulaPrior()].
#' @param updateData Data frame with node columns for the update cohort.
#' @param nRep Number of simulation replicates (default 500).
#' @param statistic Test statistic: `"pairwiseTau"` (default).
#' @param seed Optional random seed for reproducibility.
#' @return List with `obsStat`, `pValue`, `nullDistribution`, and `nRep`.
#' @export
TestPriorConsistency <- function(priorFit,
                               updateData,
                               nRep = 500L,
                               statistic = c("pairwiseTau"),
                               seed = NULL) {
  statistic <- match.arg(statistic)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  node_cols <- priorFit$nodeCols
  implied_cor <- priorFit$impliedCor
  n_update <- nrow(updateData)

  uniform_obs <- ApplyMarginalSpec(updateData, priorFit$marginalSpec, nodeCols = node_cols)
  obs_stat <- pairwise_tau_stat(uniform_obs, implied_cor)

  null_dist <- vapply(seq_len(nRep), function(i) {
    sim_scores <- SimulateFactorCopula(priorFit, nSim = n_update)
    uniform_sim <- UniformMarginalTransform(sim_scores)
    pairwise_tau_stat(uniform_sim, implied_cor)
  }, numeric(1))

  p_value <- mean(null_dist >= obs_stat)

  list(
    obsStat = obs_stat,
    pValue = p_value,
    nullDistribution = null_dist,
    nRep = nRep,
    statistic = statistic
  )
}
