# Fit Phase 2 copula update on a small sample

Applies prior marginal transforms, then fits a graphical model and/or
vine copula on the update cohort.

## Usage

``` r
FitCopulaUpdate(
  updateData,
  priorFit,
  method = c("graphical", "vine", "both"),
  nlambda = 40,
  glassoMethod = c("stars", "ebic"),
  starsThresh = 0.1
)
```

## Arguments

- updateData:

  Data frame with node columns for the update cohort.

- priorFit:

  Output of
  [`FitFactorCopulaPrior()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitFactorCopulaPrior.md).

- method:

  One of `"graphical"`, `"vine"`, or `"both"`.

- nlambda:

  Number of lambda values for glasso (graphical path).

- glassoMethod:

  Lambda selection for graphical path.

- starsThresh:

  StARS threshold.

## Value

List with `graphical` and/or `vine` sub-results, plus `uniformMatrix`.
