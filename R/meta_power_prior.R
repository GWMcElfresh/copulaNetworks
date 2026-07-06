.require_powerprior <- function() {
  if (!requireNamespace("powerprior", quietly = TRUE)) {
    stop(
      "Package 'powerprior' is required. Install with install.packages('powerprior').",
      call. = FALSE
    )
  }
}

.cohort_to_normal_scores <- function(cohort_data, marginal_spec, node_cols) {
  uniform_matrix <- ApplyMarginalSpec(cohort_data, marginal_spec, nodeCols = node_cols)
  z_matrix <- qnorm(uniform_matrix)
  colnames(z_matrix) <- node_cols
  z_matrix[, node_cols, drop = FALSE]
}

.sequential_power_prior <- function(z_list, a0) {
  .require_powerprior()
  if (length(z_list) < 1L) {
    stop("At least one cohort is required.", call. = FALSE)
  }
  if (length(a0) == 1L) {
    a0 <- rep(a0, length(z_list))
  }
  if (length(a0) != length(z_list)) {
    stop("a0 length must be 1 or match the number of cohorts.", call. = FALSE)
  }
  # ponytail: sequential NIW chaining, not joint multi-study product likelihood
  pp <- powerprior::powerprior_multivariate(z_list[[1]], a0 = a0[1])
  if (length(z_list) > 1L) {
    for (k in 2:length(z_list)) {
      pp <- powerprior::powerprior_multivariate(
        z_list[[k]],
        a0 = a0[k],
        mu0 = pp$mu_n,
        kappa0 = pp$kappa_n,
        nu0 = pp$nu_n,
        Lambda0 = pp$Lambda_n
      )
    }
  }
  pp
}

.cor_from_niw <- function(niw_obj) {
  lambda_mat <- if (!is.null(niw_obj$Lambda_star)) {
    niw_obj$Lambda_star
  } else {
    niw_obj$Lambda_n
  }
  nu_val <- if (!is.null(niw_obj$nu_star)) {
    niw_obj$nu_star
  } else {
    niw_obj$nu_n
  }
  p <- ncol(lambda_mat)
  denom <- max(nu_val - p - 1, 1e-6)
  # ponytail: NIW mean covariance -> cor; upgrade path = sample posterior cor elements
  sigma_hat <- lambda_mat / denom
  ridge <- max(1e-4, 1e-6 * mean(diag(sigma_hat), na.rm = TRUE))
  sigma_adj <- sigma_hat + diag(ridge, p)
  suppressWarnings(stats::cov2cor(sigma_adj))
}

#' Fit a meta-analytic correlation prior from multiple historical cohorts
#'
#' Pools large historical cohorts via sequential multivariate Normal-Inverse-Wishart
#' power priors ([powerprior::powerprior_multivariate()]). Requires optional
#' **powerprior** package.
#'
#' @param historicalCohorts Named list of data frames, each with `nodeCols`.
#' @param nodeCols Character vector of variable names shared across cohorts.
#' @param a0 Per-cohort discount in `[0, 1]`; recycled scalar default `0.5`.
#' @return List with `impliedCor`, `marginalSpec`, `nodeCols`, `powerPrior`,
#'   `cohortSummaries`, `zCohortList`, and fields compatible with [FitCopulaUpdate()].
#' @export
FitMetaCorPrior <- function(historicalCohorts, nodeCols, a0 = NULL) {
  .require_powerprior()
  if (!is.list(historicalCohorts) || is.null(names(historicalCohorts))) {
    stop("historicalCohorts must be a named list of data frames.", call. = FALSE)
  }
  if (length(historicalCohorts) < 1L) {
    stop("historicalCohorts must contain at least one cohort.", call. = FALSE)
  }

  n_cohorts <- length(historicalCohorts)
  if (is.null(a0)) {
    a0 <- rep(0.5, n_cohorts)
  }
  if (length(a0) == 1L) {
    a0 <- rep(a0, n_cohorts)
  }
  if (length(a0) != n_cohorts) {
    stop("a0 length must be 1 or match the number of cohorts.", call. = FALSE)
  }
  if (any(a0 < 0 | a0 > 1)) {
    stop("Each a0 must lie in [0, 1].", call. = FALSE)
  }

  pooled_data <- do.call(rbind, historicalCohorts)
  marginal_spec <- FitMarginalSpec(pooled_data, nodeCols = nodeCols)

  z_list <- lapply(historicalCohorts, function(cohort_df) {
    cohort_mat <- as.matrix(cohort_df[, nodeCols, drop = FALSE])
    if (nrow(cohort_mat) < 2L) {
      stop("Each historical cohort needs at least 2 rows for multivariate power prior.",
           call. = FALSE)
    }
    .cohort_to_normal_scores(cohort_df, marginal_spec, nodeCols)
  })
  names(z_list) <- names(historicalCohorts)

  pp <- .sequential_power_prior(z_list, a0)
  implied_cor <- .cor_from_niw(pp)
  colnames(implied_cor) <- rownames(implied_cor) <- nodeCols

  cohort_summaries <- data.frame(
    cohort = names(historicalCohorts),
    n = vapply(historicalCohorts, nrow, integer(1)),
    a0 = a0,
    stringsAsFactors = FALSE
  )

  list(
    impliedCor = implied_cor,
    marginalSpec = marginal_spec,
    nodeCols = nodeCols,
    n = sum(cohort_summaries$n),
    nHistorical = pp$kappa_n,
    powerPrior = pp,
    cohortSummaries = cohort_summaries,
    zCohortList = z_list,
    a0 = a0,
    loadings = NULL,
    nFactors = NULL
  )
}

#' Bayesian update of a meta-analytic prior with new cohort data
#'
#' Applies [powerprior::posterior_multivariate()] on normal scores from the
#' update cohort mapped through the pooled marginal specification.
#'
#' @param metaPrior Output of [FitMetaCorPrior()].
#' @param updateData Data frame with node columns for the update cohort.
#' @return List with `posterior`, `impliedCor`, and `zMatrix`.
#' @export
FitBayesianMetaUpdate <- function(metaPrior, updateData) {
  .require_powerprior()
  node_cols <- metaPrior$nodeCols
  z_update <- .cohort_to_normal_scores(updateData, metaPrior$marginalSpec, node_cols)
  posterior <- powerprior::posterior_multivariate(metaPrior$powerPrior, z_update)
  implied_cor <- .cor_from_niw(posterior)
  colnames(implied_cor) <- rownames(implied_cor) <- node_cols
  list(
    posterior = posterior,
    impliedCor = implied_cor,
    zMatrix = z_update
  )
}

#' Sensitivity of posterior correlation to global discount parameter
#'
#' Re-fits the sequential power prior at each scalar `a0` (recycled across cohorts)
#' and updates with the new cohort. ponytail: global discount only in this grid.
#'
#' @param metaPrior Output of [FitMetaCorPrior()].
#' @param updateData Update cohort data frame.
#' @param a0Values Numeric vector of discount values in `[0, 1]`.
#' @return List with `summary` (data.frame) and `corMatrices` (named list).
#' @export
RunDiscountSensitivity <- function(metaPrior, updateData, a0Values = c(0, 0.25, 0.5, 0.75, 1)) {
  .require_powerprior()
  if (is.null(metaPrior$zCohortList)) {
    stop("metaPrior$zCohortList is required (re-run FitMetaCorPrior()).", call. = FALSE)
  }
  z_list <- metaPrior$zCohortList
  node_cols <- metaPrior$nodeCols
  n_cohorts <- length(z_list)

  cor_matrices <- vector("list", length(a0Values))
  names(cor_matrices) <- as.character(a0Values)
  summary_rows <- vector("list", length(a0Values))

  for (i in seq_along(a0Values)) {
    a0_val <- a0Values[i]
    a0_vec <- rep(a0_val, n_cohorts)
    pp <- .sequential_power_prior(z_list, a0_vec)
    temp_prior <- list(
      powerPrior = pp,
      marginalSpec = metaPrior$marginalSpec,
      nodeCols = node_cols
    )
    bayes_fit <- FitBayesianMetaUpdate(temp_prior, updateData)
    cor_matrices[[i]] <- bayes_fit$impliedCor
    summary_rows[[i]] <- data.frame(
      a0 = a0_val,
      nHistorical = pp$kappa_n,
      posteriorKappa = bayes_fit$posterior$kappa_star,
      stringsAsFactors = FALSE
    )
  }
  names(cor_matrices) <- as.character(a0Values)

  list(
    summary = do.call(rbind, summary_rows),
    corMatrices = cor_matrices
  )
}

#' Run the full Bayesian meta-analysis power-prior pipeline
#'
#' Phase 1 pools historical cohorts via power prior; Phase 2 updates with
#' graphical copula fit and conjugate Bayesian posterior.
#'
#' @param historicalCohorts Named list of large reference data frames.
#' @param updateData Small update data frame.
#' @param nodeCols Character vector of node column names.
#' @param a0 Per-cohort discount parameters (passed to [FitMetaCorPrior()]).
#' @param a0Sensitivity Optional vector of global `a0` values for sensitivity analysis.
#' @param phase2Method One of `"graphical"`, `"vine"`, or `"both"`.
#' @param nlambda Glasso path length for graphical update.
#' @param glassoMethod Lambda selection for graphical update.
#' @param starsThresh StARS threshold.
#' @param outDir Optional directory for RDS checkpoints.
#' @return List with `metaPrior`, `updateFit`, `bayesFit`, `comparison`,
#'   optional `sensitivity`, and `meta`.
#' @export
RunMetaAnalysisPipeline <- function(historicalCohorts,
                                    updateData,
                                    nodeCols,
                                    a0 = NULL,
                                    a0Sensitivity = NULL,
                                    phase2Method = c("graphical", "vine", "both"),
                                    nlambda = 40,
                                    glassoMethod = c("stars", "ebic"),
                                    starsThresh = 0.1,
                                    outDir = NULL) {
  phase2_method <- match.arg(phase2Method)
  glasso_method <- match.arg(glassoMethod)
  t_start <- Sys.time()

  meta_prior <- FitMetaCorPrior(historicalCohorts, nodeCols = nodeCols, a0 = a0)
  update_fit <- FitCopulaUpdate(
    updateData,
    priorFit = meta_prior,
    method = phase2_method,
    nlambda = nlambda,
    glassoMethod = glasso_method,
    starsThresh = starsThresh
  )
  bayes_fit <- FitBayesianMetaUpdate(meta_prior, updateData)

  comparison <- NULL
  if (!is.null(update_fit$graphical)) {
    comparison <- ComparePriorToUpdate(meta_prior, update_fit)
  }

  sensitivity <- NULL
  if (!is.null(a0Sensitivity)) {
    sensitivity <- RunDiscountSensitivity(meta_prior, updateData, a0Values = a0Sensitivity)
  }

  meta <- list(
    elapsedSec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
    phase2Method = phase2_method,
    packageVersion = utils::packageVersion("copulaNetworks")
  )

  result <- list(
    metaPrior = meta_prior,
    updateFit = update_fit,
    bayesFit = bayes_fit,
    comparison = comparison,
    sensitivity = sensitivity,
    meta = meta
  )

  if (!is.null(outDir)) {
    SaveCheckpoint(result, outDir, filename = "meta_analysis_result.rds")
  }

  result
}
