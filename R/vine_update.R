#' Group variables by factor loading structure
#'
#' Heuristic clustering by loading sign and magnitude. Full factor-tree vine
#' truncation is deferred (ponytail ceiling).
#'
#' @param loadings Named numeric vector of factor loadings (Kendall taus).
#' @param k Number of groups (default: number of distinct sign buckets, max 4).
#' @return Named integer vector of group assignments.
#' @export
FactorGroupsFromLoadings <- function(loadings, k = NULL) {
  if (is.null(k)) {
    k <- min(4L, max(1L, length(unique(sign(loadings)))))
  }
  loading_df <- data.frame(
    name = names(loadings),
    value = as.numeric(loadings),
    stringsAsFactors = FALSE
  )
  if (k == 1L || nrow(loading_df) <= k) {
    groups <- setNames(rep(1L, length(loadings)), names(loadings))
    return(groups)
  }
  km <- stats::kmeans(loading_df$value, centers = k, nstart = 10)
  groups <- setNames(as.integer(km$cluster), names(loadings))
  groups
}

#' Fit a vine copula update on uniform pseudo-observations (Phase 2a)
#'
#' Requires the optional **VineCopula** package.
#'
#' @param uniformMatrix Matrix of pseudo-observations in (0, 1) (n x d).
#' @param factorGroups Optional named integer vector of factor groups (metadata).
#' @param selectionCrit Selection criterion for pair-copula families.
#' @return List with `vineFit`, `logLik`, `n`, `factorGroups`, and `familySet`.
#' @export
FitVineCopulaUpdate <- function(uniformMatrix,
                                factorGroups = NULL,
                                selectionCrit = c("AIC", "BIC")) {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop(
      "Package 'VineCopula' is required for FitVineCopulaUpdate(). ",
      "Install with install.packages('VineCopula').",
      call. = FALSE
    )
  }
  selection_crit <- match.arg(selectionCrit)
  u_matrix <- as.matrix(uniformMatrix)
  n_obs <- nrow(u_matrix)

  vine_fit <- VineCopula::RVineStructureSelect(
    data = u_matrix,
    familyset = NA,
    selectioncrit = selection_crit,
    progress = FALSE
  )

  list(
    vineFit = vine_fit,
    logLik = vine_fit$loglik,
    n = n_obs,
    factorGroups = factorGroups,
    selectionCrit = selection_crit
  )
}
