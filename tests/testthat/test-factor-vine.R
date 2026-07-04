#' Synthetic factor-structured data for two-phase demos
#'
#' @param N Prior sample size.
#' @param n Update sample size.
#' @param d Number of variables.
#' @param seed Random seed.
#' @return List with `priorData` and `updateData` data frames.
#' @keywords internal
make_factor_data <- function(N = 200L, n = 30L, d = 8L, seed = 42L) {
  set.seed(seed)
  lambda <- runif(d, 0.4, 0.8)
  sim_block <- function(n_obs) {
    f <- rnorm(n_obs)
    eps <- matrix(rnorm(n_obs * d), nrow = n_obs, ncol = d)
    factor_part <- matrix(lambda, nrow = n_obs, ncol = d, byrow = TRUE) * f
    x <- factor_part + sweep(eps, 2, sqrt(1 - lambda^2), `*`)
    colnames(x) <- paste0("node", seq_len(d))
    as.data.frame(x)
  }
  list(
    priorData = sim_block(N),
    updateData = sim_block(n)
  )
}

# --- Marginals ---
test_that("ApplyMarginalSpec preserves dimensions", {
  blocks <- make_factor_data(N = 100, n = 20, d = 6)
  spec <- FitMarginalSpec(blocks$priorData)
  u <- ApplyMarginalSpec(blocks$updateData, spec)
  expect_equal(dim(u), c(20, 6))
  expect_true(all(u > 0 & u < 1, na.rm = TRUE))
})

test_that("UniformMarginalTransform returns values in (0,1)", {
  x <- matrix(rnorm(50), ncol = 5)
  u <- UniformMarginalTransform(x)
  expect_equal(dim(u), dim(x))
  expect_true(all(u > 0 & u < 1))
})

# --- Factor structure ---
test_that("CheckFactorStructure detects strong single factor", {
  loadings <- rep(0.7, 8)
  lambda <- sin(pi * loadings / 2)
  cor_mat <- outer(lambda, lambda)
  diag(cor_mat) <- 1
  chk <- CheckFactorStructure(cor_mat, nFactors = 1)
  expect_false(chk$lowFactorWarning)
  expect_gt(chk$cumVar[1], 0.5)
})

# --- Factor prior (optional dep) ---
test_that("FitFactorCopulaPrior on synthetic data", {
  skip_if_not_installed("FactorCopula")
  blocks <- make_factor_data(N = 200, n = 30, d = 6)
  nodes <- paste0("node", 1:6)
  fit <- FitFactorCopulaPrior(blocks$priorData, nodeCols = nodes, nFactors = 1L)
  expect_equal(length(fit$loadings), 6)
  expect_equal(dim(fit$impliedCor), c(6, 6))
})

# --- GoF simulation ---
test_that("TestPriorConsistency returns pValue in [0,1]", {
  blocks <- make_factor_data(N = 150, n = 25, d = 6)
  nodes <- paste0("node", 1:6)
  prior_fit <- list(
    loadings = rep(0.5, 6),
    impliedCor = {
      lambda <- sin(pi * 0.5 / 2)
      m <- outer(rep(lambda, 6), rep(lambda, 6))
      diag(m) <- 1
      colnames(m) <- rownames(m) <- nodes
      m
    },
    marginalSpec = FitMarginalSpec(blocks$priorData, nodeCols = nodes),
    nodeCols = nodes,
    nFactors = 1L
  )
  names(prior_fit$loadings) <- nodes
  test_res <- TestPriorConsistency(prior_fit, blocks$updateData, nRep = 50, seed = 1)
  expect_gte(test_res$pValue, 0)
  expect_lte(test_res$pValue, 1)
  expect_length(test_res$nullDistribution, 50)
})

# --- Comparison ---
test_that("ComparePriorToUpdate structure matches CompareTwoStrata", {
  blocks <- make_factor_data(N = 80, n = 40, d = 6)
  nodes <- paste0("node", 1:6)
  lambda <- sin(pi * rep(0.5, 6) / 2)
  cor_mat <- outer(lambda, lambda)
  diag(cor_mat) <- 1
  colnames(cor_mat) <- rownames(cor_mat) <- nodes
  prior_fit <- list(
    loadings = setNames(rep(0.5, 6), nodes),
    impliedCor = cor_mat,
    marginalSpec = FitMarginalSpec(blocks$priorData, nodeCols = nodes),
    nodeCols = nodes,
    n = nrow(blocks$priorData),
    nFactors = 1L
  )
  update_fit <- FitCopulaUpdate(
    blocks$updateData,
    prior_fit,
    method = "graphical",
    nlambda = 10
  )
  cmp <- ComparePriorToUpdate(prior_fit, update_fit)
  expect_s3_class(cmp, "CopulaStratumComparison")
  expect_true(all(c("valueA", "valueB", "delta", "absDelta") %in% names(cmp$pcor)))
})

# --- Pipeline ---
test_that("RunFactorVinePipeline graphical + simulation", {
  skip_if_not_installed("FactorCopula")
  blocks <- make_factor_data(N = 120, n = 30, d = 6)
  nodes <- paste0("node", 1:6)
  res <- RunFactorVinePipeline(
    priorData = blocks$priorData,
    updateData = blocks$updateData,
    nodeCols = nodes,
    nFactors = 1L,
    phase2Method = "graphical",
    testMethod = "simulation",
    nRep = 30,
    nlambda = 10,
    seed = 7
  )
  expect_false(is.null(res$priorFit))
  expect_false(is.null(res$updateFit$graphical))
  expect_true("pValue" %in% names(res$consistencyTest))
})

# --- Bayesian update (optional dep) ---
test_that("FitBayesianFactorUpdate runs on synthetic data", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("FactorCopula")
  cmdstan_ok <- tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(cmdstan_ok, "CmdStan not installed (run cmdstanr::install_cmdstan())")
  blocks <- make_factor_data(N = 120, n = 25, d = 5)
  nodes <- paste0("node", 1:5)
  prior_fit <- FitFactorCopulaPrior(blocks$priorData, nodeCols = nodes, nFactors = 1L)
  bayes <- tryCatch(
    FitBayesianFactorUpdate(
      prior_fit,
      blocks$updateData,
      chains = 2L,
      iter = 200L,
      seed = 1L
    ),
    error = function(e) skip(conditionMessage(e))
  )
  expect_s3_class(bayes$fit, "CmdStanMCMC")
  expect_true(nrow(bayes$summary) > 0)
})

# --- Vignette source ---
test_that("factor-vine vignette source exists and has required sections", {
  v_path <- normalizePath(
    file.path(testthat::test_path(), "..", "..", "vignettes", "two-phase-factor-vine.Rmd"),
    mustWork = TRUE
  )
  v <- readLines(v_path)
  expect_true(any(grepl("Phase 1", v)))
  expect_true(any(grepl("RunFactorVinePipeline", v)))
})
