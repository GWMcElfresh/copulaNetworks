# Fit marginal specification from a prior (large) sample

Stores empirical rank knots per column so update data can be mapped
through the prior ECDF. Covariate-adjusted marginals are deferred -
supply pre-adjusted residuals in `priorMatrix`.

## Usage

``` r
FitMarginalSpec(priorMatrix, nodeCols = colnames(priorMatrix))
```

## Arguments

- priorMatrix:

  Numeric matrix or data frame (N x d) from the prior cohort.

- nodeCols:

  Character vector of column names (default: all columns).

## Value

List with `method`, `nRef`, `knots`, and `nodeCols`.
