# Build an approximate correlation matrix from factor loadings

Converts Kendall tau loadings to Pearson correlations via the Gaussian
copula identity, then assembles a one- or two-factor correlation
structure.

## Usage

``` r
ImpliedFactorCorrelation(loadings, nFactors = 1L)
```

## Arguments

- loadings:

  Named numeric vector of Kendall tau loadings per variable.

- nFactors:

  Number of factors used in fitting (1 or 2).

## Value

Symmetric correlation matrix.
