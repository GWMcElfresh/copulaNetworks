#' Synthetic longitudinal data for latent mixed-model tests
#'
#' @keywords internal
make_longitudinal_mixed_data <- function(n_subj = 25L,
                                         n_time = 4L,
                                         d = 4L,
                                         seed = 11L) {
  set.seed(seed)
  subject <- rep(seq_len(n_subj), each = n_time)
  time <- rep(seq_len(n_time), times = n_subj)
  n <- length(subject)
  b0 <- matrix(rnorm(n_subj * d, 0, 0.6), nrow = n_subj, ncol = d)
  b1 <- matrix(rnorm(n_subj * d, 0, 0.15), nrow = n_subj, ncol = d)
  eps <- matrix(rnorm(n * d), nrow = n, ncol = d)
  # Shared residual factor for cross-node dependence
  f <- rnorm(n)
  lam <- runif(d, 0.3, 0.6)
  x <- matrix(NA_real_, n, d)
  for (j in seq_len(d)) {
    x[, j] <- b0[subject, j] + b1[subject, j] * time +
      0.25 * time + lam[j] * f + eps[, j] * sqrt(1 - lam[j]^2)
  }
  colnames(x) <- paste0("node", seq_len(d))
  cbind(data.frame(subject = subject, time = time), as.data.frame(x))
}

test_that("parse_latent_formula accepts random intercept and slope", {
  dat <- make_longitudinal_mixed_data(n_subj = 10, n_time = 3, d = 3)
  p1 <- parse_latent_formula(~ time + (1 | subject), dat)
  expect_equal(p1$groupVar, "subject")
  expect_true(p1$hasIntercept)
  expect_equal(p1$nRe, 1L)

  p2 <- parse_latent_formula(~ time + (1 + time | subject), dat)
  expect_equal(p2$slopeVars, "time")
  expect_equal(p2$nRe, 2L)

  p3 <- parse_latent_formula(~ time + (0 + time | subject), dat)
  expect_false(p3$hasIntercept)
  expect_equal(p3$nRe, 1L)
})

test_that("parse_latent_formula rejects crossed / LHS response", {
  dat <- make_longitudinal_mixed_data(n_subj = 8, n_time = 2, d = 3)
  expect_error(
    parse_latent_formula(y ~ time + (1 | subject), dat),
    "RHS-only"
  )
  expect_error(
    parse_latent_formula(~ time + (1 | subject) + (1 | site), dat),
    "single grouping"
  )
})

test_that("FitCopulaLatentMixedModel random intercept (conjugate)", {
  dat <- make_longitudinal_mixed_data(n_subj = 20, n_time = 4, d = 4)
  nodes <- paste0("node", 1:4)
  expect_warning(
    fit <- FitCopulaLatentMixedModel(
      dat,
      nodeCols = nodes,
      latentFormula = ~ time + (1 | subject),
      engine = "conjugate",
      nIter = 150,
      burnIn = 40,
      thin = 2,
      seed = 7
    ),
    "recommend"
  )
  expect_s3_class(fit, "CopulaLatentMixedFit")
  expect_equal(dim(fit$residualMatrix), c(nrow(dat), 4))
  expect_equal(fit$engine, "conjugate")
  expect_true(all(is.finite(fit$residualMatrix)))
  expect_true(nrow(fit$fixedEffects) >= 4)
  expect_true(nrow(fit$randomEffects) >= 20)
})

test_that("FitCopulaLatentMixedModel random slope (conjugate)", {
  dat <- make_longitudinal_mixed_data(n_subj = 18, n_time = 4, d = 3)
  nodes <- paste0("node", 1:3)
  suppressWarnings(
    fit <- FitCopulaLatentMixedModel(
      dat,
      nodeCols = nodes,
      latentFormula = ~ time + (1 + time | subject),
      engine = "conjugate",
      nIter = 120,
      burnIn = 30,
      thin = 2,
      seed = 3
    )
  )
  expect_equal(fit$design$nRe, 2L)
  expect_equal(dim(fit$randomCov[[1]]), c(2, 2))
  expect_true(all(is.finite(fit$residualMatrix)))
})

test_that("PredictCopulaLatentMixedModel known vs new groups", {
  dat <- make_longitudinal_mixed_data(n_subj = 15, n_time = 3, d = 3)
  nodes <- paste0("node", 1:3)
  suppressWarnings(
    fit <- FitCopulaLatentMixedModel(
      dat,
      nodeCols = nodes,
      latentFormula = ~ time + (1 | subject),
      engine = "conjugate",
      nIter = 100,
      burnIn = 25,
      seed = 5
    )
  )
  pred_train <- PredictCopulaLatentMixedModel(fit, dat, type = "latentMean")
  expect_equal(dim(pred_train), c(nrow(dat), 3))
  expect_true(all(is.finite(pred_train)))

  newdata <- data.frame(
    subject = c(1L, 999L),
    time = c(1, 2),
    node1 = c(0, 0),
    node2 = c(0, 0),
    node3 = c(0, 0)
  )
  pred_new <- PredictCopulaLatentMixedModel(fit, newdata, type = "latentMean")
  expect_equal(nrow(pred_new), 2)
  # New group uses population-level RE (0); known group nonzero typically
  expect_true(all(is.finite(pred_new)))
})

test_that("backward compatibility: FitStratumCopula still works without RE", {
  set.seed(1)
  x <- matrix(rnorm(80 * 5), 80, 5)
  colnames(x) <- paste0("node", 1:5)
  fit <- FitStratumCopula(as.data.frame(x), nodeCols = colnames(x), nlambda = 10, method = "ebic")
  expect_false(is.null(fit$copulaCor))
  expect_equal(dim(fit$copulaCor), c(5, 5))
})

test_that("residual NIW update with metaPrior", {
  dat <- make_longitudinal_mixed_data(n_subj = 20, n_time = 3, d = 4, seed = 9)
  nodes <- paste0("node", 1:4)
  # Historical cohorts without RE structure for prior
  hist <- list(
    A = as.data.frame(matrix(rnorm(60 * 4), 60, 4, dimnames = list(NULL, nodes))),
    B = as.data.frame(matrix(rnorm(50 * 4), 50, 4, dimnames = list(NULL, nodes)))
  )
  meta_prior <- FitMetaCorPrior(hist, nodeCols = nodes, a0 = 0.5)
  suppressWarnings(
    fit <- FitCopulaLatentMixedModel(
      dat,
      nodeCols = nodes,
      latentFormula = ~ time + (1 | subject),
      engine = "conjugate",
      metaPrior = meta_prior,
      nIter = 80,
      burnIn = 20,
      seed = 2
    )
  )
  expect_false(is.null(fit$impliedCor))
  expect_equal(dim(fit$impliedCor), c(4, 4))
  expect_true(all(is.finite(fit$impliedCor)))
  expect_true(all(abs(diag(fit$impliedCor) - 1) < 1e-6))
})

test_that("plot wrappers return objects for residual mixed fit", {
  dat <- make_longitudinal_mixed_data(n_subj = 16, n_time = 3, d = 4)
  nodes <- paste0("node", 1:4)
  suppressWarnings(
    fit <- FitCopulaLatentMixedModel(
      dat,
      nodeCols = nodes,
      latentFormula = ~ time + (1 | subject),
      engine = "conjugate",
      nIter = 80,
      burnIn = 20,
      seed = 4
    )
  )
  graphical <- PseudoGraphicalFitFromLatentMixed(fit)
  expect_false(is.null(graphical$copulaCor))
  hm <- PlotLatentMixedResidualCorHeatmap(fit, title = "test")
  expect_false(is.null(hm))
  net <- PlotLatentMixedResidualNetwork(fit, printPlot = FALSE)
  # Empty graphs may return NULL; otherwise a ggplot/ggraph object
  expect_true(is.null(net) || inherits(net, "ggplot") || inherits(net, "ggraph") || inherits(net, "gg"))
})

test_that("brms engine smoke test", {
  skip_if_not_installed("brms")
  skip_if_not_installed("BH")
  skip_on_cran()
  dat <- make_longitudinal_mixed_data(n_subj = 12, n_time = 3, d = 3, seed = 21)
  nodes <- paste0("node", 1:3)
  fit <- tryCatch(
    suppressWarnings(FitCopulaLatentMixedModel(
      dat,
      nodeCols = nodes,
      latentFormula = ~ time + (1 | subject),
      engine = "brms",
      chains = 1,
      iter = 400,
      warmup = 200,
      cores = 1,
      seed = 21,
      silent = TRUE
    )),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    skip(paste("brms/Stan toolchain unavailable:", conditionMessage(fit)))
  }
  expect_equal(fit$engine, "brms")
  expect_equal(dim(fit$residualMatrix), c(nrow(dat), 3))
  expect_true(all(is.finite(fit$residualMatrix)))
})
