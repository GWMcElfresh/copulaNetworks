#' Draw from Inverse-Wishart (scale parameterization matching NIW literature)
#'
#' @param nu Degrees of freedom.
#' @param Psi Scale matrix (posterior scale).
#' @return Draw of Sigma ~ IW(nu, Psi).
#' @keywords internal
#' @noRd
draw_inv_wishart <- function(nu, Psi) {
  p <- nrow(Psi)
  Psi_inv <- tryCatch(
    solve(Psi),
    error = function(e) solve(Psi + diag(1e-8, p))
  )
  # W ~ Wishart(nu, Psi^{-1}) => W^{-1} ~ IW(nu, Psi)
  W <- stats::rWishart(1L, df = max(nu, p), Sigma = Psi_inv)[, , 1L]
  solve(W)
}

#' Conjugate Gibbs sampler for one node's latent mixed model
#'
#' Model: z = X beta + Z_re b_g + eps, eps ~ N(0, sigma2),
#' b_g ~ N(0, D). Flat prior on beta; IG on sigma2; IW/IG on D.
#'
#' @param z Numeric response (latent scores for one node).
#' @param design Output of [build_latent_design_matrices()].
#' @param nIter MCMC iterations.
#' @param burnIn Burn-in iterations.
#' @param thin Thinning interval.
#' @param seed Optional RNG seed.
#' @return List with posterior means for beta, b, D, sigma2, residuals, samples.
#' @keywords internal
#' @noRd
fit_latent_mixed_conjugate_node <- function(z,
                                            design,
                                            nIter = 800L,
                                            burnIn = 200L,
                                            thin = 2L,
                                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  X <- design$X
  Z <- design$Z
  group_id <- design$groupId
  n_obs <- length(z)
  p_fixed <- ncol(X)
  q_re <- design$nRe
  n_groups <- design$nGroups

  # Priors
  a0_sig <- 0.01
  b0_sig <- 0.01
  nu0_d <- q_re + 2
  Psi0_d <- diag(q_re)

  beta <- stats::coef(stats::lm.fit(X, z))
  beta[!is.finite(beta)] <- 0
  b_mat <- matrix(0, nrow = n_groups, ncol = q_re)
  colnames(b_mat) <- design$reNames
  rownames(b_mat) <- design$groupLevels
  D <- diag(q_re)
  sigma2 <- stats::var(z)
  if (!is.finite(sigma2) || sigma2 < 1e-8) {
    sigma2 <- 1
  }

  keep_idx <- seq(burnIn + 1L, nIter, by = thin)
  n_keep <- length(keep_idx)
  beta_samps <- matrix(NA_real_, nrow = n_keep, ncol = p_fixed)
  colnames(beta_samps) <- design$fixedNames
  D_samps <- array(NA_real_, dim = c(q_re, q_re, n_keep))
  sigma2_samps <- numeric(n_keep)
  b_sum <- matrix(0, nrow = n_groups, ncol = q_re)
  resid_sum <- numeric(n_obs)
  store_i <- 0L

  XtX <- crossprod(X)
  for (iter in seq_len(nIter)) {
    # Linear predictor from RE
    re_contrib <- expand_random_effects(b_mat, group_id, Z)
    y_tilde <- z - re_contrib

    # beta | rest ~ N(m, sigma2 * (X'X)^{-1}) with flat prior
    XtX_inv <- tryCatch(
      solve(XtX),
      error = function(e) solve(XtX + diag(1e-8, p_fixed))
    )
    beta_hat <- as.numeric(XtX_inv %*% crossprod(X, y_tilde))
    beta <- as.numeric(
      beta_hat + sqrt(sigma2) * (t(chol(XtX_inv)) %*% stats::rnorm(p_fixed))
    )

    fixed_contrib <- as.numeric(X %*% beta)
    y_re <- z - fixed_contrib

    # b_g | rest
    D_inv <- tryCatch(
      solve(D),
      error = function(e) solve(D + diag(1e-8, q_re))
    )
    for (g in seq_len(n_groups)) {
      idx <- which(group_id == g)
      Zg <- Z[idx, , drop = FALSE]
      yg <- y_re[idx]
      prec <- D_inv + crossprod(Zg) / sigma2
      V <- tryCatch(solve(prec), error = function(e) solve(prec + diag(1e-8, q_re)))
      m <- as.numeric(V %*% (crossprod(Zg, yg) / sigma2))
      b_mat[g, ] <- as.numeric(m + t(chol(V)) %*% stats::rnorm(q_re))
    }

    re_contrib <- expand_random_effects(b_mat, group_id, Z)
    resid <- z - fixed_contrib - re_contrib

    # sigma2 | rest ~ IG
    sse <- sum(resid^2)
    a_n <- a0_sig + n_obs / 2
    b_n <- b0_sig + sse / 2
    sigma2 <- 1 / stats::rgamma(1, shape = a_n, rate = b_n)

    # D | rest ~ IW / IG
    S_b <- crossprod(b_mat)
    if (q_re == 1L) {
      # Inverse-Gamma as IW(nu, Psi) with p=1
      nu_n <- nu0_d + n_groups
      psi_n <- as.numeric(Psi0_d) + as.numeric(S_b)
      D <- matrix(1 / stats::rgamma(1, shape = nu_n / 2, rate = psi_n / 2), 1, 1)
    } else {
      nu_n <- nu0_d + n_groups
      Psi_n <- Psi0_d + S_b
      D <- draw_inv_wishart(nu_n, Psi_n)
    }
    colnames(D) <- rownames(D) <- design$reNames

    if (iter %in% keep_idx) {
      store_i <- store_i + 1L
      beta_samps[store_i, ] <- beta
      D_samps[, , store_i] <- D
      sigma2_samps[store_i] <- sigma2
      b_sum <- b_sum + b_mat
      resid_sum <- resid_sum + resid
    }
  }

  beta_mean <- colMeans(beta_samps)
  b_mean <- b_sum / n_keep
  colnames(b_mean) <- design$reNames
  rownames(b_mean) <- design$groupLevels
  D_mean <- apply(D_samps, c(1, 2), mean)
  colnames(D_mean) <- rownames(D_mean) <- design$reNames
  resid_mean <- resid_sum / n_keep

  list(
    beta = beta_mean,
    randomEffects = b_mean,
    randomCov = D_mean,
    sigma2 = mean(sigma2_samps),
    residuals = resid_mean,
    engine = "conjugate",
    samples = list(
      beta = beta_samps,
      randomCov = D_samps,
      sigma2 = sigma2_samps
    )
  )
}

#' Fit all nodes with the conjugate Gibbs engine
#'
#' @param zMatrix n x p latent score matrix.
#' @param design Shared design matrices.
#' @param ... Passed to [fit_latent_mixed_conjugate_node()].
#' @return Named list of per-node fits.
#' @keywords internal
#' @noRd
fit_latent_mixed_conjugate <- function(zMatrix, design, ...) {
  warning(
    "Using engine = 'conjugate' (Gibbs). ",
    "We recommend engine = 'brms' (Hamiltonian Monte Carlo) when brms is available.",
    call. = FALSE
  )
  node_cols <- colnames(zMatrix)
  fits <- lapply(node_cols, function(nm) {
    fit_latent_mixed_conjugate_node(zMatrix[, nm], design = design, ...)
  })
  names(fits) <- node_cols
  fits
}
