# Group variables by factor loading structure

Heuristic clustering by loading sign and magnitude. Full factor-tree
vine truncation is deferred (ponytail ceiling).

## Usage

``` r
FactorGroupsFromLoadings(loadings, k = NULL)
```

## Arguments

- loadings:

  Named numeric vector of factor loadings (Kendall taus).

- k:

  Number of groups (default: number of distinct sign buckets, max 4).

## Value

Named integer vector of group assignments.
