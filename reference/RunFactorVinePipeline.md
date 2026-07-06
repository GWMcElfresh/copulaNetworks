# Run the full two-phase factor-vine copula pipeline

Phase 1 fits a factor copula prior on a large reference cohort; Phase 2
updates with graphical and/or vine models on a small sample; optional
consistency testing via simulation and/or Bayesian updating.

## Usage

``` r
RunFactorVinePipeline(
  priorData,
  updateData,
  nodeCols,
  nFactors = 1L,
  linkingCopula = "bvn",
  nQuad = 25L,
  phase2Method = c("graphical", "vine", "both"),
  testMethod = c("simulation", "bayes", "both", "none"),
  nRep = 500L,
  nlambda = 40,
  glassoMethod = c("stars", "ebic"),
  starsThresh = 0.1,
  outDir = NULL,
  seed = NULL
)
```

## Arguments

- priorData:

  Large reference data frame (N x variables).

- updateData:

  Small update data frame (n x variables).

- nodeCols:

  Character vector of node column names.

- nFactors:

  Number of latent factors for Phase 1 (1 or 2).

- linkingCopula:

  Linking copula family per variable (Phase 1).

- nQuad:

  Quadrature points for FactorCopula integration.

- phase2Method:

  One of `"graphical"`, `"vine"`, or `"both"`.

- testMethod:

  One of `"simulation"`, `"bayes"`, `"both"`, or `"none"`.

- nRep:

  Simulation replicates for
  [`TestPriorConsistency()`](https://GWMcElfresh.github.io/copulaNetworks/reference/TestPriorConsistency.md).

- nlambda:

  Glasso path length for graphical update.

- glassoMethod:

  Lambda selection for graphical update.

- starsThresh:

  StARS threshold.

- outDir:

  Optional directory for RDS checkpoints.

- seed:

  Optional random seed for simulation test.

## Value

List with `priorFit`, `updateFit`, `comparison`, `consistencyTest`, and
`meta`.
