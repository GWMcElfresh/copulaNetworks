# Compare pairwise matrix values across two strata

Compare pairwise matrix values across two strata

## Usage

``` r
ComparePairwiseMatrices(
  resA,
  resB,
  labelA,
  labelB,
  matrixName = c("pcor", "copulaCor"),
  edgeOnly = NULL,
  deltaThreshold = 0.05
)
```

## Arguments

- resA:

  First stratum fit result.

- resB:

  Second stratum fit result.

- labelA:

  Label for stratum A.

- labelB:

  Label for stratum B.

- matrixName:

  One of `"pcor"` or `"copulaCor"`.

- edgeOnly:

  For `"pcor"`, if `TRUE` only adjacency edges are non-zero. For
  `"copulaCor"`, defaults to `FALSE` (all pairs).

- deltaThreshold:

  Threshold for direction classification.

## Value

Data frame with comparison columns.
