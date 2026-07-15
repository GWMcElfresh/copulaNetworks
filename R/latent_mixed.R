#' Map observed node columns to latent Gaussian scores
#'
#' @param data Data frame.
#' @param nodeCols Node column names.
#' @param marginalSpec Optional prior marginal specification.
#' @return n x p matrix of normal scores.
#' @keywords internal
#' @noRd
latent_scores_from_data <- function(data, nodeCols, marginalSpec = NULL) {
  missing_cols <- setdiff(nodeCols, names(data))
  if (length(missing_cols)) {
    stop("Missing nodeCols in data: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.null(marginalSpec)) {
    uniform_matrix <- ApplyMarginalSpec(data, marginalSpec, nodeCols = nodeCols)
    z_matrix <- qnorm(uniform_matrix)
  } else {
    input_matrix <- as.matrix(data[, nodeCols, drop = FALSE])
    z_matrix <- NonparanormalTransform(input_matrix)
  }
  colnames(z_matrix) <- nodeCols
  z_matrix
}

#' Assemble camelCase summaries from per-node fits
#'
#' @keywords internal
#' @noRd
assemble_latent_mixed_summaries <- function(node_fits, design) {
  node_cols <- names(node_fits)
  fixed_effects <- do.call(rbind, lapply(node_cols, function(nm) {
    beta <- node_fits[[nm]]$beta
    data.frame(
      node = nm,
      term = names(beta),
      estimate = as.numeric(beta),
      stringsAsFactors = FALSE
    )
  }))

  random_effects <- do.call(rbind, lapply(node_cols, function(nm) {
    b_mat <- node_fits[[nm]]$randomEffects
    out <- as.data.frame(as.table(b_mat))
    names(out) <- c("group", "term", "estimate")
    out$node <- nm
    out$groupVar <- design$groupVar
    out[, c("node", "groupVar", "group", "term", "estimate")]
  }))

  random_cov <- lapply(node_fits, function(f) f$randomCov)
  residual_matrix <- do.call(cbind, lapply(node_fits, function(f) f$residuals))
  colnames(residual_matrix) <- node_cols

  list(
    fixedEffects = fixed_effects,
    randomEffects = random_effects,
    randomCov = random_cov,
    residualMatrix = residual_matrix
  )
}

#' Fit Bayesian copula latent mixed models
#'
#' Random effects are introduced **only** in the latent Gaussian layer after
#' nonparanormal (or prior-ECDF) transformation of `nodeCols`. A shared
#' RHS-only `latentFormula` specifies the mean structure applied to each node's
#' latent scores \eqn{z_j}:
#'
#' \deqn{z_j = X\beta_j + Z b_j + \varepsilon_j, \quad
#' b_j \sim \mathrm{MVN}(0, D_j), \quad
#' \varepsilon \sim \mathrm{MVN}(0, \Sigma)}
#'
#' Cross-node dependence is captured by residual \eqn{\Sigma}. When `metaPrior`
#' from [FitMetaCorPrior()] is supplied, residual scores are passed to the
#' existing NIW / powerprior update ([FitBayesianMetaUpdate()]) without changing
#' the conjugate NIW machinery.
#'
#' Supported random-effect structures (single grouping factor): `(1 | g)`,
#' `(0 + x | g)`, `(1 + x | g)` / `(x | g)`.
#'
#' @param data Data frame with node columns and covariates / grouping variables.
#' @param nodeCols Character vector of variables entering the copula.
#' @param latentFormula RHS-only formula for the **latent** mean structure,
#'   e.g. `~ time + treatment + (1 + time | subject)`. Does not select which
#'   variables enter the copula (`nodeCols` does).
#' @param engine `"brms"` (recommended HMC) or `"conjugate"` (Gibbs; always
#'   available). Conjugate emits a warning recommending brms.
#' @param metaPrior Optional output of [FitMetaCorPrior()] for residual NIW update.
#' @param marginalSpec Optional [FitMarginalSpec()] result; if `NULL` and
#'   `metaPrior` is supplied, uses `metaPrior$marginalSpec`; else same-cohort NPN.
#' @param nIter,burnIn,thin Conjugate MCMC controls.
#' @param chains,iter,warmup,cores brms MCMC controls (ignored for conjugate).
#' @param seed Optional RNG seed.
#' @param silent If `TRUE`, suppress brms progress output.
#' @return List with camelCase fields: `nodeCols`, `latentFormula`, `engine`,
#'   `zMatrix`, `residualMatrix`, `nodeFits`, `fixedEffects`, `randomEffects`,
#'   `randomCov`, `design`, `marginalSpec`, and optionally `bayesFit` /
#'   `impliedCor` when `metaPrior` is supplied. Also `residualGraphical` when
#'   a glasso fit on residuals succeeds.
#' @export
#' @examples
#' set.seed(1)
#' n_subj <- 20
#' n_time <- 4
#' dat <- data.frame(
#'   subject = rep(seq_len(n_subj), each = n_time),
#'   time = rep(seq_len(n_time), times = n_subj)
#' )
#' b0 <- rnorm(n_subj, 0, 0.5)
#' dat$y1 <- b0[dat$subject] + 0.3 * dat$time + rnorm(nrow(dat))
#' dat$y2 <- 0.7 * dat$y1 + rnorm(nrow(dat), sd = 0.5)
#' dat$y3 <- 0.4 * dat$y1 + rnorm(nrow(dat), sd = 0.7)
#' fit <- FitCopulaLatentMixedModel(
#'   dat,
#'   nodeCols = c("y1", "y2", "y3"),
#'   latentFormula = ~ time + (1 | subject),
#'   engine = "conjugate",
#'   nIter = 200,
#'   burnIn = 50
#' )
#' dim(fit$residualMatrix)
FitCopulaLatentMixedModel <- function(data,
                                      nodeCols,
                                      latentFormula,
                                      engine = c("brms", "conjugate"),
                                      metaPrior = NULL,
                                      marginalSpec = NULL,
                                      nIter = 800L,
                                      burnIn = 200L,
                                      thin = 2L,
                                      chains = 2L,
                                      iter = 1000L,
                                      warmup = 500L,
                                      cores = 1L,
                                      seed = NULL,
                                      silent = TRUE) {
  engine <- match.arg(engine)
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  if (!length(nodeCols) || !is.character(nodeCols)) {
    stop("nodeCols must be a non-empty character vector.", call. = FALSE)
  }

  if (is.null(marginalSpec) && !is.null(metaPrior)) {
    marginalSpec <- metaPrior$marginalSpec
  }

  parsed <- parse_latent_formula(latentFormula, data)
  design <- build_latent_design_matrices(parsed, data)
  z_matrix <- latent_scores_from_data(data, nodeCols, marginalSpec = marginalSpec)

  if (engine == "brms") {
    node_fits <- fit_latent_mixed_brms(
      zMatrix = z_matrix,
      data = data,
      latentFormula = latentFormula,
      design = design,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      silent = silent
    )
  } else {
    node_fits <- fit_latent_mixed_conjugate(
      zMatrix = z_matrix,
      design = design,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      seed = seed
    )
  }

  summaries <- assemble_latent_mixed_summaries(node_fits, design)

  bayes_fit <- NULL
  implied_cor <- NULL
  if (!is.null(metaPrior)) {
    if (is.null(metaPrior$powerPrior)) {
      stop("metaPrior must contain powerPrior (from FitMetaCorPrior()).", call. = FALSE)
    }
    # Residuals are already on the latent scale; call powerprior NIW directly.
    posterior <- powerprior::posterior_multivariate(
      metaPrior$powerPrior,
      as.matrix(summaries$residualMatrix)
    )
    implied_cor <- .cor_from_niw(
      posterior,
      fallback = stats::cor(summaries$residualMatrix)
    )
    colnames(implied_cor) <- rownames(implied_cor) <- nodeCols
    bayes_fit <- list(
      posterior = posterior,
      impliedCor = implied_cor,
      zMatrix = summaries$residualMatrix
    )
  }

  residual_graphical <- tryCatch(
    {
      resid_df <- as.data.frame(summaries$residualMatrix)
      FitStratumCopula(
        resid_df,
        nodeCols = nodeCols,
        nlambda = 20,
        method = "ebic",
        preTransformed = TRUE
      )
    },
    error = function(e) NULL
  )

  result <- list(
    nodeCols = nodeCols,
    latentFormula = latentFormula,
    engine = engine,
    zMatrix = z_matrix,
    residualMatrix = summaries$residualMatrix,
    nodeFits = node_fits,
    fixedEffects = summaries$fixedEffects,
    randomEffects = summaries$randomEffects,
    randomCov = summaries$randomCov,
    design = design,
    marginalSpec = marginalSpec,
    metaPrior = metaPrior,
    bayesFit = bayes_fit,
    impliedCor = implied_cor,
    residualGraphical = residual_graphical
  )
  class(result) <- c("CopulaLatentMixedFit", class(result))
  result
}

#' Predict latent means or residuals from a copula latent mixed fit
#'
#' @param fit Output of [FitCopulaLatentMixedModel()].
#' @param newdata Data frame with covariates / grouping variable (and node
#'   columns when `type = "residual"`).
#' @param type `"latentMean"` for \eqn{X\beta_j + Z b_j}; `"residual"` for
#'   \eqn{z_j} minus latent mean when node columns are present.
#' @return Matrix (n x p) of predictions named by `fit$nodeCols`.
#' @export
PredictCopulaLatentMixedModel <- function(fit,
                                          newdata,
                                          type = c("latentMean", "residual")) {
  type <- match.arg(type)
  if (!inherits(fit, "CopulaLatentMixedFit") && is.null(fit$nodeFits)) {
    stop("fit must be output from FitCopulaLatentMixedModel().", call. = FALSE)
  }
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame.", call. = FALSE)
  }

  parsed <- parse_latent_formula(fit$latentFormula, newdata)
  design_new <- build_latent_design_matrices(parsed, newdata)
  node_cols <- fit$nodeCols
  n_obs <- nrow(newdata)
  pred <- matrix(NA_real_, nrow = n_obs, ncol = length(node_cols))
  colnames(pred) <- node_cols

  for (nm in node_cols) {
    node_fit <- fit$nodeFits[[nm]]
    beta <- node_fit$beta
    # Align beta to new X columns
    beta_vec <- stats::setNames(rep(0, ncol(design_new$X)), colnames(design_new$X))
    shared_beta <- intersect(names(beta), names(beta_vec))
    beta_vec[shared_beta] <- beta[shared_beta]
    fixed_contrib <- as.numeric(design_new$X %*% beta_vec)

    b_mat <- node_fit$randomEffects
    b_new <- matrix(
      0,
      nrow = design_new$nGroups,
      ncol = fit$design$nRe,
      dimnames = list(design_new$groupLevels, fit$design$reNames)
    )
    known <- intersect(rownames(b_mat), rownames(b_new))
    shared_re <- intersect(colnames(b_mat), colnames(b_new))
    if (length(known) && length(shared_re)) {
      b_new[known, shared_re] <- b_mat[known, shared_re]
    }
    # New groups remain 0 (population level)
    re_contrib <- expand_random_effects(b_new, design_new$groupId, design_new$Z)
    pred[, nm] <- fixed_contrib + re_contrib
  }

  if (type == "residual") {
    missing_nodes <- setdiff(node_cols, names(newdata))
    if (length(missing_nodes)) {
      stop(
        "type = 'residual' requires node columns in newdata. Missing: ",
        paste(missing_nodes, collapse = ", "),
        call. = FALSE
      )
    }
    z_new <- latent_scores_from_data(newdata, node_cols, marginalSpec = fit$marginalSpec)
    pred <- z_new - pred
  }

  pred
}
