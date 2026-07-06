# Map update values through a prior marginal specification

Map update values through a prior marginal specification

## Usage

``` r
ApplyMarginalSpec(updateMatrix, marginalSpec, nodeCols = marginalSpec$nodeCols)
```

## Arguments

- updateMatrix:

  Numeric matrix or data frame (n x d) from the update cohort.

- marginalSpec:

  Output of
  [`FitMarginalSpec()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitMarginalSpec.md).

- nodeCols:

  Character vector of columns to transform.

## Value

Matrix of pseudo-observations in (0, 1).
