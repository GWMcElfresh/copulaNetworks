#' Nonparanormal transformation (rank-based Gaussian copula marginals)
#'
#' Maps each column to normal scores via ranks:
#' \deqn{\hat{Z}_j = \Phi^{-1}(\mathrm{rank}(X_j) / (n + 1))}
#'
#' @param inputMatrix Numeric matrix (n x p).
#' @return Matrix of normal scores with column names preserved.
#' @export
NonparanormalTransform <- function(inputMatrix) {
  n <- nrow(inputMatrix)
  quantile_matrix <- apply(inputMatrix, 2, function(x) {
    qnorm(rank_to_uniform(x, n_ref = n))
  })
  colnames(quantile_matrix) <- colnames(inputMatrix)
  quantile_matrix
}

#' Fit a nonparanormal graphical model on one stratum
#'
#' @param data Data frame containing node columns.
#' @param nodeCols Character vector of endogenous variable names.
#' @param nlambda Number of lambda values for the glasso path.
#' @param method Lambda selection criterion: `"stars"` (StARS) or `"ebic"`.
#' @param starsThresh StARS stability threshold (used when `method = "stars"`).
#' @param preTransformed If `TRUE`, `data` columns are already normal scores
#'   (skip [NonparanormalTransform()]). Ignored when `quantileMatrix` is supplied.
#' @param quantileMatrix Optional pre-computed normal-score matrix (n x p).
#' @return List with correlation matrix, partial correlation matrix, selected
#'   graph, etc. Returns `NULL` if fewer than 3 non-constant variables.
#' @export
FitStratumCopula <- function(data,
                             nodeCols,
                             nlambda = 40,
                             method = c("stars", "ebic"),
                             starsThresh = 0.1,
                             preTransformed = FALSE,
                             quantileMatrix = NULL) {
  method <- match.arg(method)
  input_matrix <- as.matrix(data[, nodeCols, drop = FALSE])

  column_variances <- apply(input_matrix, 2, var, na.rm = TRUE)
  keep <- column_variances > 1e-10
  if (sum(keep) < 3) {
    warning("Fewer than 3 non-constant variables in stratum — skipping.", call. = FALSE)
    return(NULL)
  }
  input_matrix <- input_matrix[, keep, drop = FALSE]
  kept_cols <- colnames(input_matrix)

  quantile_matrix <- if (!is.null(quantileMatrix)) {
    quantileMatrix[, kept_cols, drop = FALSE]
  } else if (isTRUE(preTransformed)) {
    input_matrix
  } else {
    NonparanormalTransform(input_matrix)
  }
  copula_cor <- cor(quantile_matrix, use = "pairwise.complete.obs")

  fit <- huge::huge(quantile_matrix, method = "glasso", nlambda = nlambda, verbose = FALSE)

  if (method == "stars") {
    sel <- huge::huge.select(fit, criterion = "stars", stars.thresh = starsThresh, verbose = FALSE)
  } else {
    sel <- huge::huge.select(fit, criterion = "ebic", verbose = FALSE)
  }

  Omega <- as.matrix(sel$opt.icov)
  colnames(Omega) <- rownames(Omega) <- kept_cols

  D_inv <- diag(1 / sqrt(diag(Omega)))
  pcor <- -D_inv %*% Omega %*% D_inv
  diag(pcor) <- 1
  colnames(pcor) <- rownames(pcor) <- kept_cols

  adjacency_matrix <- (abs(pcor) > 1e-8) * 1
  diag(adjacency_matrix) <- 0

  list(
    copulaCor = copula_cor,
    precision = Omega,
    pcor = pcor,
    adjacency = adjacency_matrix,
    keptCols = kept_cols,
    n = nrow(input_matrix),
    lambdaOpt = sel$opt.lambda,
    quantileMatrix = quantile_matrix
  )
}

#' Fit copula models across prepared strata
#'
#' @param prep Output of [PrepareCopulaData()].
#' @param method Lambda selection criterion.
#' @param nlambda Number of lambda values.
#' @param starsThresh StARS threshold.
#' @param minN Minimum observations required per stratum.
#' @param includeFull If `TRUE`, also fit on the full (unstratified) dataset.
#' @param outDir Optional directory to save `fits.rds`.
#' @return Named list of fit results keyed by stratum name.
#' @export
FitCopulaStrata <- function(prep,
                            method = c("stars", "ebic"),
                            nlambda = 40,
                            starsThresh = 0.1,
                            minN = 10,
                            includeFull = FALSE,
                            outDir = NULL) {
  method <- match.arg(method)
  node_cols <- prep$nodeCols
  fits <- list()

  for (nm in names(prep$strata)) {
    stratum_data <- prep$strata[[nm]]
    n_obs <- nrow(stratum_data)
    if (n_obs < minN) {
      warning("Skipping stratum '", nm, "' — only ", n_obs, " observations (minN = ", minN, ").",
              call. = FALSE)
      next
    }
    message("Fitting copula for stratum: ", nm, " (n = ", n_obs, ")")
    stratum_result <- FitStratumCopula(
      stratum_data,
      nodeCols = node_cols,
      nlambda = nlambda,
      method = method,
      starsThresh = starsThresh
    )
    if (!is.null(stratum_result)) {
      fits[[nm]] <- stratum_result
    }
  }

  if (isTRUE(includeFull)) {
    message("Fitting copula for full dataset (n = ", nrow(prep$data), ")")
    res_full <- FitStratumCopula(
      prep$data,
      nodeCols = node_cols,
      nlambda = nlambda,
      method = method,
      starsThresh = starsThresh
    )
    if (!is.null(res_full)) {
      fits[["__full__"]] <- res_full
    }
  }

  result <- list(
    fits = fits,
    nodeCols = node_cols,
    method = method,
    nlambda = nlambda,
    starsThresh = starsThresh,
    minN = minN
  )

  if (!is.null(outDir)) {
    SaveCheckpoint(result, outDir, filename = "fits.rds")
  }

  result
}
