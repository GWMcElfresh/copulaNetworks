#' Nonparanormal transformation (rank-based Gaussian copula marginals)
#'
#' Maps each column to normal scores via ranks:
#' \deqn{\hat{Z}_j = \Phi^{-1}(\mathrm{rank}(X_j) / (n + 1))}
#'
#' @param X Numeric matrix (n x p).
#' @return Matrix of normal scores with column names preserved.
#' @export
nonparanormal_transform <- function(X) {
  n <- nrow(X)
  Z <- apply(X, 2, function(x) {
    r <- rank(x, ties.method = "average")
    qnorm(r / (n + 1))
  })
  colnames(Z) <- colnames(X)
  Z
}

#' Fit a nonparanormal graphical model on one stratum
#'
#' @param data Data frame containing node columns.
#' @param node_cols Character vector of endogenous variable names.
#' @param nlambda Number of lambda values for the glasso path.
#' @param method Lambda selection criterion: `"stars"` (StARS) or `"ebic"`.
#' @param stars_thresh StARS stability threshold (used when `method = "stars"`).
#' @return List with correlation matrix, partial correlation matrix, selected
#'   graph, etc. Returns `NULL` if fewer than 3 non-constant variables.
#' @export
fit_stratum_copula <- function(data,
                               node_cols,
                               nlambda = 40,
                               method = c("stars", "ebic"),
                               stars_thresh = 0.1) {
  method <- match.arg(method)
  X <- as.matrix(data[, node_cols, drop = FALSE])

  col_vars <- apply(X, 2, var, na.rm = TRUE)
  keep <- col_vars > 1e-10
  if (sum(keep) < 3) {
    warning("Fewer than 3 non-constant variables in stratum — skipping.", call. = FALSE)
    return(NULL)
  }
  X <- X[, keep, drop = FALSE]
  kept_cols <- colnames(X)

  Z <- nonparanormal_transform(X)
  copula_cor <- cor(Z, use = "pairwise.complete.obs")

  fit <- huge::huge(Z, method = "glasso", nlambda = nlambda, verbose = FALSE)

  if (method == "stars") {
    sel <- huge::huge.select(fit, criterion = "stars", stars.thresh = stars_thresh, verbose = FALSE)
  } else {
    sel <- huge::huge.select(fit, criterion = "ebic", verbose = FALSE)
  }

  Omega <- as.matrix(sel$opt.icov)
  colnames(Omega) <- rownames(Omega) <- kept_cols

  D_inv <- diag(1 / sqrt(diag(Omega)))
  pcor <- -D_inv %*% Omega %*% D_inv
  diag(pcor) <- 1
  colnames(pcor) <- rownames(pcor) <- kept_cols

  adj <- (abs(pcor) > 1e-8) * 1
  diag(adj) <- 0

  list(
    copula_cor = copula_cor,
    precision = Omega,
    pcor = pcor,
    adjacency = adj,
    kept_cols = kept_cols,
    n = nrow(X),
    lambda_opt = sel$opt.lambda,
    Z = Z
  )
}

#' Fit copula models across prepared strata
#'
#' @param prep Output of [prepare_copula_data()].
#' @param method Lambda selection criterion.
#' @param nlambda Number of lambda values.
#' @param stars_thresh StARS threshold.
#' @param min_n Minimum observations required per stratum.
#' @param include_full If `TRUE`, also fit on the full (unstratified) dataset.
#' @param out_dir Optional directory to save `fits.rds`.
#' @return Named list of fit results keyed by stratum name.
#' @export
fit_copula_strata <- function(prep,
                              method = c("stars", "ebic"),
                              nlambda = 40,
                              stars_thresh = 0.1,
                              min_n = 10,
                              include_full = FALSE,
                              out_dir = NULL) {
  method <- match.arg(method)
  node_cols <- prep$node_cols
  fits <- list()

  for (nm in names(prep$strata)) {
    stratum_data <- prep$strata[[nm]]
    n_obs <- nrow(stratum_data)
    if (n_obs < min_n) {
      warning("Skipping stratum '", nm, "' — only ", n_obs, " observations (min_n = ", min_n, ").",
              call. = FALSE)
      next
    }
    message("Fitting copula for stratum: ", nm, " (n = ", n_obs, ")")
    res <- fit_stratum_copula(
      stratum_data,
      node_cols = node_cols,
      nlambda = nlambda,
      method = method,
      stars_thresh = stars_thresh
    )
    if (!is.null(res)) {
      fits[[nm]] <- res
    }
  }

  if (isTRUE(include_full)) {
    message("Fitting copula for full dataset (n = ", nrow(prep$data), ")")
    res_full <- fit_stratum_copula(
      prep$data,
      node_cols = node_cols,
      nlambda = nlambda,
      method = method,
      stars_thresh = stars_thresh
    )
    if (!is.null(res_full)) {
      fits[["__full__"]] <- res_full
    }
  }

  result <- list(
    fits = fits,
    node_cols = node_cols,
    method = method,
    nlambda = nlambda,
    stars_thresh = stars_thresh,
    min_n = min_n
  )

  if (!is.null(out_dir)) {
    save_checkpoint(result, out_dir, filename = "fits.rds")
  }

  result
}
