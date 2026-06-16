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

test_that("nonparanormal_transform preserves dimensions", {
  X <- matrix(rnorm(50), ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  Z <- nonparanormal_transform(X)
  expect_equal(dim(Z), dim(X))
  expect_equal(colnames(Z), colnames(X))
})

test_that("fit_stratum_copula returns symmetric pcor with unit diagonal", {
  dat <- make_synthetic_data()
  res <- fit_stratum_copula(dat, node_cols = paste0("node", 1:8), nlambda = 10)
  expect_false(is.null(res))
  expect_equal(dim(res$pcor), c(length(res$kept_cols), length(res$kept_cols)))
  expect_equal(unname(diag(res$pcor)), rep(1, ncol(res$pcor)), tolerance = 1e-6)
  expect_equal(res$pcor, t(res$pcor), tolerance = 1e-6)
})

test_that("build_strata respects filter and min_n", {
  dat <- make_synthetic_data()
  strata <- build_strata(
    dat,
    spec = list(
      filter = quote(stratum == "A"),
      group_by = "stratum",
      min_n = 5
    )
  )
  expect_equal(length(strata), 1)
  expect_equal(nrow(strata[[1]]), 30)
})

test_that("prepare_copula_data builds strata from specs", {
  dat <- make_synthetic_data()
  prep <- prepare_copula_data(
    data = dat,
    id_cols = "SubjectId",
    strata_cols = c("stratum", "Vaccine"),
    node_cols = paste0("node", 1:8),
    strata_specs = list(
      by_stratum = list(stratum_col = "stratum")
    )
  )
  expect_equal(length(prep$strata), 2)
  expect_equal(length(prep$node_cols), 8)
})

test_that("fit_copula_strata fits all strata", {
  dat <- make_synthetic_data()
  prep <- prepare_copula_data(
    data = dat,
    id_cols = "SubjectId",
    strata_cols = c("stratum", "Vaccine"),
    node_cols = paste0("node", 1:8),
    strata_specs = list(by_stratum = list(stratum_col = "stratum"))
  )
  fits <- fit_copula_strata(prep, nlambda = 10, min_n = 10)
  expect_equal(length(fits$fits), 2)
})

test_that("compare_two_strata returns both matrix types", {
  dat <- make_synthetic_data()
  prep <- prepare_copula_data(
    data = dat,
    id_cols = "SubjectId",
    strata_cols = c("stratum", "Vaccine"),
    node_cols = paste0("node", 1:8),
    strata_specs = list(by_stratum = list(stratum_col = "stratum"))
  )
  fits <- fit_copula_strata(prep, nlambda = 10, min_n = 10)
  cmp <- compare_two_strata(
    fits$fits[["by_stratum::A"]],
    fits$fits[["by_stratum::B"]],
    label_a = "A",
    label_b = "B"
  )
  expect_true("pcor" %in% names(cmp))
  expect_true("copula_cor" %in% names(cmp))
  expect_true(all(c("value_a", "value_b", "delta", "abs_delta") %in% names(cmp$pcor)))
})

test_that("checkpoint save and load roundtrip", {
  tmp <- tempfile("copula_ckpt")
  dir.create(tmp)
  obj <- list(x = 1:3)
  save_checkpoint(obj, tmp, filename = "test.rds")
  loaded <- load_checkpoint(tmp, filename = "test.rds")
  expect_equal(loaded, obj)
})
