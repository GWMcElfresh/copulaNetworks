#' Synthetic factor-structured cohorts for meta-analysis demos
#'
#' @param d Number of variables.
#' @param cohort_sizes Named vector of cohort sample sizes.
#' @param n_update Update cohort size.
#' @param seed Random seed.
#' @return List with `historicalCohorts` and `updateData`.
#' @keywords internal
make_meta_cohorts <- function(d = 6L,
                              cohort_sizes = c(A = 80L, B = 70L),
                              n_update = 25L,
                              seed = 42L) {
  set.seed(seed)
  lambda <- runif(d, 0.4, 0.75)
  sim_block <- function(n_obs, jitter = 0) {
    lam <- pmax(0.2, pmin(0.85, lambda + jitter))
    f <- rnorm(n_obs)
    eps <- matrix(rnorm(n_obs * d), nrow = n_obs, ncol = d)
    factor_part <- matrix(lam, nrow = n_obs, ncol = d, byrow = TRUE) * f
    x <- factor_part + sweep(eps, 2, sqrt(1 - lam^2), `*`)
    colnames(x) <- paste0("node", seq_len(d))
    as.data.frame(x)
  }
  historical <- lapply(seq_along(cohort_sizes), function(i) {
    sim_block(cohort_sizes[i], jitter = 0.02 * (i - 1))
  })
  names(historical) <- names(cohort_sizes)
  list(
    historicalCohorts = historical,
    updateData = sim_block(n_update)
  )
}

test_that("FitMetaCorPrior returns expected structure", {
  blocks <- make_meta_cohorts(d = 6)
  nodes <- paste0("node", 1:6)
  fit <- FitMetaCorPrior(blocks$historicalCohorts, nodeCols = nodes, a0 = c(0.5, 0.5))
  expect_equal(dim(fit$impliedCor), c(6, 6))
  expect_equal(fit$nodeCols, nodes)
  expect_equal(nrow(fit$cohortSummaries), 2)
  expect_false(is.null(fit$zCohortList))
  expect_true(all(abs(diag(fit$impliedCor) - 1) < 1e-8))
})

test_that("FitBayesianMetaUpdate returns posterior correlation", {
  blocks <- make_meta_cohorts(d = 6)
  nodes <- paste0("node", 1:6)
  meta_prior <- FitMetaCorPrior(blocks$historicalCohorts, nodeCols = nodes)
  bayes <- FitBayesianMetaUpdate(meta_prior, blocks$updateData)
  expect_equal(dim(bayes$impliedCor), c(6, 6))
  expect_false(is.null(bayes$posterior$kappa_star))
})

test_that("RunDiscountSensitivity grid has expected dimensions", {
  blocks <- make_meta_cohorts(d = 6)
  nodes <- paste0("node", 1:6)
  meta_prior <- FitMetaCorPrior(blocks$historicalCohorts, nodeCols = nodes)
  a0_grid <- c(0, 0.5, 1)
  sens <- RunDiscountSensitivity(
    meta_prior,
    blocks$updateData,
    a0Values = a0_grid
  )
  expect_equal(nrow(sens$summary), length(a0_grid))
  expect_equal(length(sens$corMatrices), length(a0_grid))
  expect_equal(dim(sens$corMatrices[[1]]), c(6, 6))
})

test_that("RunDiscountSensitivity cor matrices are finite at a0 = 0", {
  blocks <- make_meta_cohorts(d = 6)
  nodes <- paste0("node", 1:6)
  meta_prior <- FitMetaCorPrior(blocks$historicalCohorts, nodeCols = nodes)
  sens <- RunDiscountSensitivity(
    meta_prior,
    blocks$updateData,
    a0Values = c(0, 0.5, 1)
  )
  expect_true(all(is.finite(sens$corMatrices[["0"]])))
})

test_that("RunMetaAnalysisPipeline runs end-to-end", {
  blocks <- make_meta_cohorts(d = 6)
  nodes <- paste0("node", 1:6)
  res <- RunMetaAnalysisPipeline(
    historicalCohorts = blocks$historicalCohorts,
    updateData = blocks$updateData,
    nodeCols = nodes,
    phase2Method = "graphical",
    nlambda = 10,
    a0Sensitivity = c(0, 1)
  )
  expect_false(is.null(res$metaPrior$impliedCor))
  expect_false(is.null(res$updateFit$graphical))
  expect_false(is.null(res$bayesFit$impliedCor))
  expect_equal(nrow(res$sensitivity$summary), 2)
})

test_that("FitCopulaUpdate works with meta prior without loadings", {
  blocks <- make_meta_cohorts(d = 6)
  nodes <- paste0("node", 1:6)
  meta_prior <- FitMetaCorPrior(blocks$historicalCohorts, nodeCols = nodes)
  update_fit <- FitCopulaUpdate(
    blocks$updateData,
    priorFit = meta_prior,
    method = "graphical",
    nlambda = 10
  )
  expect_false(is.null(update_fit$graphical))
  expect_null(update_fit$vine)
})

test_that("meta-analysis vignette HTML includes embedded figures", {
  v_path <- system.file("doc", "meta-analysis-power-prior.html", package = "copulaNetworks")
  if (!nzchar(v_path)) {
    pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/")
    candidates <- c(
      file.path(pkg_root, "inst", "doc", "meta-analysis-power-prior.html"),
      file.path(pkg_root, "vignettes", "meta-analysis-power-prior.html")
    )
    v_path <- candidates[file.exists(candidates)][1]
  }
  skip_if_not(
    length(v_path) == 1 && nzchar(v_path) && file.exists(v_path),
    "Vignette HTML not built (run tools::buildVignettes() or R CMD build)"
  )
  v <- readLines(v_path)
  expect_true(any(grepl("RunMetaAnalysisPipeline", v)))
  expect_true(
    any(grepl("<img", v)),
    info = "Vignette HTML has no embedded figures"
  )
  expect_true(
    any(grepl("Meta-analytic prior network", v)),
    info = "Prior network figure missing from vignette"
  )
})
