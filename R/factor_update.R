#' Build a pseudo graphical fit from an implied correlation matrix
#'
#' Converts an implied correlation matrix (e.g. from a factor copula prior) into
#' the structure expected by [PlotCopulaNetwork()] and the heatmap helpers.
#'
#' @param impliedCor Symmetric correlation matrix.
#' @param nObs Reference sample size for metadata.
#' @return List with `copulaCor`, `pcor`, `adjacency`, `keptCols`, and related fields.
#' @export
PseudoGraphicalFitFromCor <- function(impliedCor, nObs = NA_integer_) {
  cor_mat <- as.matrix(impliedCor)
  kept_cols <- colnames(cor_mat)
  n_vars <- ncol(cor_mat)

  if (n_vars < 3) {
    stop("Implied correlation matrix must have at least 3 variables.", call. = FALSE)
  }

  precision <- tryCatch(
    solve(cor_mat),
    error = function(e) {
      cor_adj <- cor_mat + diag(1e-4, n_vars)
      solve(cor_adj)
    }
  )
  colnames(precision) <- rownames(precision) <- kept_cols

  d_inv <- diag(1 / sqrt(diag(precision)))
  pcor <- -d_inv %*% precision %*% d_inv
  diag(pcor) <- 1
  colnames(pcor) <- rownames(pcor) <- kept_cols

  adjacency_matrix <- (abs(pcor) > 1e-8) * 1
  diag(adjacency_matrix) <- 0

  list(
    copulaCor = cor_mat,
    precision = precision,
    pcor = pcor,
    adjacency = adjacency_matrix,
    keptCols = kept_cols,
    n = nObs,
    lambdaOpt = NA_real_,
    quantileMatrix = NULL
  )
}

#' Fit Phase 2 copula update on a small sample
#'
#' Applies prior marginal transforms, then fits a graphical model and/or vine
#' copula on the update cohort.
#'
#' @param updateData Data frame with node columns for the update cohort.
#' @param priorFit Prior fit with `marginalSpec`, `impliedCor`, and `nodeCols`
#'   (from [FitFactorCopulaPrior()] or [FitMetaCorPrior()]).
#' @param method One of `"graphical"`, `"vine"`, or `"both"`.
#' @param nlambda Number of lambda values for glasso (graphical path).
#' @param glassoMethod Lambda selection for graphical path.
#' @param starsThresh StARS threshold.
#' @return List with `graphical` and/or `vine` sub-results, plus `uniformMatrix`.
#' @export
FitCopulaUpdate <- function(updateData,
                            priorFit,
                            method = c("graphical", "vine", "both"),
                            nlambda = 40,
                            glassoMethod = c("stars", "ebic"),
                            starsThresh = 0.1) {
  method <- match.arg(method)
  glasso_method <- match.arg(glassoMethod)
  node_cols <- priorFit$nodeCols
  marginal_spec <- priorFit$marginalSpec

  uniform_matrix <- ApplyMarginalSpec(updateData, marginal_spec, nodeCols = node_cols)
  z_matrix <- qnorm(uniform_matrix)

  result <- list(uniformMatrix = uniform_matrix, zMatrix = z_matrix)

  if (method %in% c("graphical", "both")) {
    update_df <- as.data.frame(z_matrix)
    result$graphical <- FitStratumCopula(
      update_df,
      nodeCols = node_cols,
      nlambda = nlambda,
      method = glasso_method,
      starsThresh = starsThresh,
      preTransformed = TRUE
    )
  }

  if (method %in% c("vine", "both")) {
    if (is.null(priorFit$loadings)) {
      if (method == "vine") {
        stop(
          "priorFit$loadings is required for vine update (e.g. FitFactorCopulaPrior()).",
          call. = FALSE
        )
      }
    } else {
      n_factors <- if (is.null(priorFit$nFactors)) 1L else priorFit$nFactors
      factor_groups <- FactorGroupsFromLoadings(priorFit$loadings, k = n_factors)
      result$vine <- FitVineCopulaUpdate(
        uniform_matrix,
        factorGroups = factor_groups
      )
    }
  }

  result
}

#' Compare prior factor model to Phase 2 update fit
#'
#' Builds a pseudo-fit from the prior implied correlation and delegates to
#' [CompareTwoStrata()].
#'
#' @param priorFit Prior fit with `impliedCor` and `n` ([FitFactorCopulaPrior()]
#'   or [FitMetaCorPrior()]).
#' @param updateFit Output of [FitCopulaUpdate()].
#' @param deltaThreshold Threshold for direction labels.
#' @return A `CopulaStratumComparison` object.
#' @export
ComparePriorToUpdate <- function(priorFit, updateFit, deltaThreshold = 0.05) {
  if (is.null(updateFit$graphical)) {
    stop("updateFit$graphical is required for ComparePriorToUpdate().", call. = FALSE)
  }
  prior_pseudo <- PseudoGraphicalFitFromCor(priorFit$impliedCor, nObs = priorFit$n)
  CompareTwoStrata(
    prior_pseudo,
    updateFit$graphical,
    labelA = "prior",
    labelB = "update",
    deltaThreshold = deltaThreshold
  )
}
