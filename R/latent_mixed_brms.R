.require_brms <- function() {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for engine = 'brms'. ",
      "Install with install.packages('brms'), or use engine = 'conjugate'.",
      call. = FALSE
    )
  }
}

#' Build a brms formula for one latent node
#'
#' @param node_name Response column name in the working data frame.
#' @param latentFormula Original RHS-only latent formula.
#' @return A brms-compatible formula.
#' @keywords internal
#' @noRd
build_brms_node_formula <- function(node_name, latentFormula) {
  rhs <- paste(deparse(latentFormula[[2]]), collapse = " ")
  stats::as.formula(paste(node_name, "~", rhs))
}

#' Fit one node with brms HMC
#'
#' @param z Numeric latent scores.
#' @param data Covariate data frame (same rows as z).
#' @param node_name Temporary response name.
#' @param latentFormula RHS-only formula.
#' @param design Design list from [build_latent_design_matrices()].
#' @param chains,iter,warmup,cores,seed brms MCMC controls.
#' @param silent Suppress brms progress.
#' @return List matching conjugate node-fit fields.
#' @keywords internal
#' @noRd
fit_latent_mixed_brms_node <- function(z,
                                       data,
                                       node_name,
                                       latentFormula,
                                       design,
                                       chains = 2L,
                                       iter = 1000L,
                                       warmup = 500L,
                                       cores = 1L,
                                       seed = NULL,
                                       silent = TRUE) {
  .require_brms()
  work <- data
  work[[node_name]] <- z
  form <- build_brms_node_formula(node_name, latentFormula)

  brms_args <- list(
    formula = form,
    data = work,
    family = stats::gaussian(),
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    refresh = if (isTRUE(silent)) 0 else NULL
  )
  # Drop NULL refresh for older brms
  if (is.null(brms_args$refresh)) {
    brms_args$refresh <- NULL
  }

  fit <- suppressWarnings(do.call(brms::brm, brms_args))

  # Fixed effects: prefer summary()$fixed (brms); fall back to coef()
  beta_aligned <- stats::setNames(rep(0, length(design$fixedNames)), design$fixedNames)
  fixed_tab <- tryCatch(as.data.frame(summary(fit)$fixed), error = function(e) NULL)
  if (!is.null(fixed_tab) && nrow(fixed_tab) > 0L) {
    est_col <- if ("Estimate" %in% names(fixed_tab)) "Estimate" else 1L
    for (nm in rownames(fixed_tab)) {
      target <- nm
      if (nm == "Intercept") target <- "(Intercept)"
      if (target %in% names(beta_aligned)) {
        beta_aligned[target] <- as.numeric(fixed_tab[nm, est_col])
      }
    }
  } else {
    beta_coef <- tryCatch(stats::coef(fit), error = function(e) NULL)
    if (!is.null(beta_coef)) {
      for (nm in names(beta_coef)) {
        target <- if (nm == "Intercept") "(Intercept)" else nm
        if (target %in% names(beta_aligned)) {
          beta_aligned[target] <- as.numeric(beta_coef[[nm]])
        }
      }
    }
  }

  ranef_list <- brms::ranef(fit)
  group_var <- design$groupVar
  if (!group_var %in% names(ranef_list)) {
    stop("brms ranef() missing grouping factor '", group_var, "'.", call. = FALSE)
  }
  re_arr <- ranef_list[[group_var]]
  # re_arr: groups x stats x effects; Estimate is typically [, , "Estimate"] or dim2
  if (length(dim(re_arr)) == 3L) {
    # groups x Estimate/Est.Error/... x effects  OR groups x effects x stats
    dn <- dimnames(re_arr)
    if (!is.null(dn[[2]]) && "Estimate" %in% dn[[2]]) {
      b_est <- re_arr[, "Estimate", , drop = FALSE]
      b_mat <- matrix(b_est, nrow = dim(re_arr)[1], ncol = dim(re_arr)[3])
      colnames(b_mat) <- dn[[3]]
      rownames(b_mat) <- dn[[1]]
    } else if (!is.null(dn[[3]]) && "Estimate" %in% dn[[3]]) {
      b_est <- re_arr[, , "Estimate", drop = FALSE]
      b_mat <- matrix(b_est, nrow = dim(re_arr)[1], ncol = dim(re_arr)[2])
      colnames(b_mat) <- dn[[2]]
      rownames(b_mat) <- dn[[1]]
    } else {
      b_mat <- as.matrix(re_arr[, , 1])
    }
  } else {
    b_mat <- as.matrix(re_arr)
  }
  # Reorder rows to design$groupLevels
  b_aligned <- matrix(
    0,
    nrow = design$nGroups,
    ncol = design$nRe,
    dimnames = list(design$groupLevels, design$reNames)
  )
  for (g in rownames(b_mat)) {
    if (g %in% rownames(b_aligned)) {
      shared <- intersect(colnames(b_mat), colnames(b_aligned))
      if (length(shared)) {
        b_aligned[g, shared] <- b_mat[g, shared]
      }
      # Map Intercept naming
      if ("Intercept" %in% colnames(b_mat) && "(Intercept)" %in% colnames(b_aligned)) {
        b_aligned[g, "(Intercept)"] <- b_mat[g, "Intercept"]
      }
    }
  }

  # Random-effect covariance from VarCorr
  vc <- brms::VarCorr(fit)
  D_mean <- diag(design$nRe)
  colnames(D_mean) <- rownames(D_mean) <- design$reNames
  if (group_var %in% names(vc) && !is.null(vc[[group_var]]$cov)) {
    cov_est <- vc[[group_var]]$cov
    if (is.list(cov_est) && !is.null(cov_est$Estimate)) {
      cov_mat <- as.matrix(cov_est$Estimate)
    } else if (is.array(cov_est) && length(dim(cov_est)) == 3L) {
      cov_mat <- cov_est[, , 1]
    } else {
      cov_mat <- as.matrix(cov_est)
    }
    shared <- intersect(colnames(cov_mat), design$reNames)
    if (length(shared)) {
      D_mean[shared, shared] <- cov_mat[shared, shared]
    }
    if ("Intercept" %in% colnames(cov_mat) && "(Intercept)" %in% design$reNames) {
      D_mean["(Intercept)", "(Intercept)"] <- cov_mat["Intercept", "Intercept"]
    }
  } else if (group_var %in% names(vc) && !is.null(vc[[group_var]]$sd)) {
    sd_est <- vc[[group_var]]$sd
    if (is.data.frame(sd_est) || is.matrix(sd_est)) {
      sds <- as.numeric(sd_est[, "Estimate"])
      nms <- rownames(sd_est)
    } else if (is.list(sd_est) && !is.null(sd_est$Estimate)) {
      sds <- as.numeric(sd_est$Estimate)
      nms <- names(sd_est$Estimate)
    } else {
      sds <- as.numeric(sd_est)
      nms <- design$reNames
    }
    if (length(sds) == 1L && design$nRe == 1L) {
      D_mean[1, 1] <- sds^2
    } else {
      for (i in seq_along(sds)) {
        nm <- if (!is.null(nms) && length(nms) >= i) nms[i] else design$reNames[i]
        if (nm == "Intercept") nm <- "(Intercept)"
        if (nm %in% design$reNames) {
          D_mean[nm, nm] <- sds[i]^2
        }
      }
    }
  }

  fixed_contrib <- as.numeric(design$X %*% beta_aligned)
  re_contrib <- expand_random_effects(b_aligned, design$groupId, design$Z)
  resid <- z - fixed_contrib - re_contrib

  sigma2 <- tryCatch(
    {
      s <- brms::VarCorr(fit)$residual$sd
      if (is.data.frame(s)) as.numeric(s["Estimate"])^2 else as.numeric(s)^2
    },
    error = function(e) stats::var(resid)
  )

  list(
    beta = beta_aligned,
    randomEffects = b_aligned,
    randomCov = D_mean,
    sigma2 = sigma2,
    residuals = resid,
    engine = "brms",
    brmsFit = fit,
    samples = NULL
  )
}

#' Fit all nodes with brms
#'
#' @keywords internal
#' @noRd
fit_latent_mixed_brms <- function(zMatrix,
                                  data,
                                  latentFormula,
                                  design,
                                  chains = 2L,
                                  iter = 1000L,
                                  warmup = 500L,
                                  cores = 1L,
                                  seed = NULL,
                                  silent = TRUE) {
  .require_brms()
  node_cols <- colnames(zMatrix)
  fits <- lapply(seq_along(node_cols), function(j) {
    nm <- node_cols[j]
    # Use a safe temporary response name
    resp <- paste0(".latent_", j)
    node_seed <- if (is.null(seed)) NULL else as.integer(seed) + j
    fit_latent_mixed_brms_node(
      z = zMatrix[, nm],
      data = data,
      node_name = resp,
      latentFormula = latentFormula,
      design = design,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = node_seed,
      silent = silent
    )
  })
  names(fits) <- node_cols
  fits
}
