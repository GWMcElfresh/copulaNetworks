#' Two-phase factor-vine copula demo (synthetic data)
#'
#' Requires optional packages FactorCopula (Phase 1) and VineCopula (vine path).
#' Run interactively: source(system.file("examples", "example_factor_vine.R",
#' package = "copulaNetworks"))

make_factor_data <- function(N = 500L, n = 30L, d = 12L, seed = 42L) {
  set.seed(seed)
  lambda <- runif(d, 0.35, 0.75)
  sim_block <- function(n_obs) {
    f <- rnorm(n_obs)
    eps <- matrix(rnorm(n_obs * d), nrow = n_obs, ncol = d)
    factor_part <- matrix(lambda, nrow = n_obs, ncol = d, byrow = TRUE) * f
    x <- factor_part + sweep(eps, 2, sqrt(1 - lambda^2), `*`)
    colnames(x) <- paste0("node", seq_len(d))
    as.data.frame(x)
  }
  list(priorData = sim_block(N), updateData = sim_block(n))
}

blocks <- make_factor_data()
nodes <- paste0("node", 1:12)

message("Phase 1: factor copula prior (requires FactorCopula)")
if (requireNamespace("FactorCopula", quietly = TRUE)) {
  prior_fit <- FitFactorCopulaPrior(
    blocks$priorData,
    nodeCols = nodes,
    nFactors = 1L,
    linkingCopula = "bvn"
  )
  print(prior_fit$loadings)

  message("Phase 2: graphical update")
  update_fit <- FitCopulaUpdate(
    blocks$updateData,
    priorFit = prior_fit,
    method = "graphical",
    nlambda = 20
  )

  cmp <- ComparePriorToUpdate(prior_fit, update_fit)
  print(head(cmp$pcor))

  message("Consistency test (simulation)")
  gof <- TestPriorConsistency(prior_fit, blocks$updateData, nRep = 100, seed = 1)
  print(gof$pValue)

  message("Full pipeline")
  res <- RunFactorVinePipeline(
    priorData = blocks$priorData,
    updateData = blocks$updateData,
    nodeCols = nodes,
    phase2Method = "graphical",
    testMethod = "simulation",
    nRep = 50,
    nlambda = 20,
    seed = 1
  )
  print(res$consistencyTest$pValue)
} else {
  message("Install FactorCopula for the full demo.")
}
