make_synthetic_data <- function(n = 60, seed = 42) {
  set.seed(seed)
  stratum <- rep(c("A", "B"), each = n / 2)
  data.frame(
    SubjectId = seq_len(n),
    stratum = stratum,
    Vaccine = ifelse(stratum == "A", "VaxA", "VaxB"),
    node1 = rnorm(n),
    node2 = rnorm(n),
    node3 = rnorm(n),
    node4 = rnorm(n),
    node5 = rnorm(n),
    node6 = rnorm(n),
    node7 = rnorm(n),
    node8 = rnorm(n),
    stringsAsFactors = FALSE
  )
}

test_that("NonparanormalTransform preserves dimensions", {
  input_matrix <- matrix(rnorm(50), ncol = 5)
  colnames(input_matrix) <- paste0("v", 1:5)
  quantile_matrix <- NonparanormalTransform(input_matrix)
  expect_equal(dim(quantile_matrix), dim(input_matrix))
  expect_equal(colnames(quantile_matrix), colnames(input_matrix))
})

test_that("FitStratumCopula returns symmetric pcor with unit diagonal", {
  dat <- make_synthetic_data()
  res <- FitStratumCopula(dat, nodeCols = paste0("node", 1:8), nlambda = 10)
  expect_false(is.null(res))
  expect_equal(dim(res$pcor), c(length(res$keptCols), length(res$keptCols)))
  expect_equal(unname(diag(res$pcor)), rep(1, ncol(res$pcor)), tolerance = 1e-6)
  expect_equal(res$pcor, t(res$pcor), tolerance = 1e-6)
})

test_that("BuildStrata respects filter and minN", {
  dat <- make_synthetic_data()
  strata <- BuildStrata(
    dat,
    spec = list(
      filter = quote(stratum == "A"),
      group_by = "stratum",
      minN = 5
    )
  )
  expect_equal(length(strata), 1)
  expect_equal(nrow(strata[[1]]), 30)
})

test_that("PrepareCopulaData builds strata from specs", {
  dat <- make_synthetic_data()
  prep <- PrepareCopulaData(
    data = dat,
    idCols = "SubjectId",
    strataCols = c("stratum", "Vaccine"),
    nodeCols = paste0("node", 1:8),
    strataSpecs = list(
      by_stratum = list(stratumCol = "stratum")
    )
  )
  expect_equal(length(prep$strata), 2)
  expect_equal(length(prep$nodeCols), 8)
})

test_that("FitCopulaStrata fits all strata", {
  dat <- make_synthetic_data()
  prep <- PrepareCopulaData(
    data = dat,
    idCols = "SubjectId",
    strataCols = c("stratum", "Vaccine"),
    nodeCols = paste0("node", 1:8),
    strataSpecs = list(by_stratum = list(stratumCol = "stratum"))
  )
  fits <- FitCopulaStrata(prep, nlambda = 10, minN = 10)
  expect_equal(length(fits$fits), 2)
})

test_that("CompareTwoStrata returns both matrix types", {
  dat <- make_synthetic_data()
  prep <- PrepareCopulaData(
    data = dat,
    idCols = "SubjectId",
    strataCols = c("stratum", "Vaccine"),
    nodeCols = paste0("node", 1:8),
    strataSpecs = list(by_stratum = list(stratumCol = "stratum"))
  )
  fits <- FitCopulaStrata(prep, nlambda = 10, minN = 10)
  cmp <- CompareTwoStrata(
    fits$fits[["by_stratum::A"]],
    fits$fits[["by_stratum::B"]],
    labelA = "A",
    labelB = "B"
  )
  expect_true("pcor" %in% names(cmp))
  expect_true("copulaCor" %in% names(cmp))
  expect_true(all(c("valueA", "valueB", "delta", "absDelta") %in% names(cmp$pcor)))
  expect_s3_class(cmp, "CopulaStratumComparison")
})

test_that("checkpoint save and load roundtrip", {
  tmp <- tempfile("copula_ckpt")
  dir.create(tmp)
  checkpoint_object <- list(x = 1:3)
  SaveCheckpoint(checkpoint_object, tmp, filename = "test.rds")
  loaded <- LoadCheckpoint(tmp, filename = "test.rds")
  expect_equal(loaded, checkpoint_object)
})
