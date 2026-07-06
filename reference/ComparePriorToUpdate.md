# Compare prior factor model to Phase 2 update fit

Builds a pseudo-fit from the prior implied correlation and delegates to
[`CompareTwoStrata()`](https://GWMcElfresh.github.io/copulaNetworks/reference/CompareTwoStrata.md).

## Usage

``` r
ComparePriorToUpdate(priorFit, updateFit, deltaThreshold = 0.05)
```

## Arguments

- priorFit:

  Output of
  [`FitFactorCopulaPrior()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitFactorCopulaPrior.md).

- updateFit:

  Output of
  [`FitCopulaUpdate()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitCopulaUpdate.md).

- deltaThreshold:

  Threshold for direction labels.

## Value

A `CopulaStratumComparison` object.
