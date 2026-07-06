# Compare two fitted strata across matrix types

Compare two fitted strata across matrix types

## Usage

``` r
CompareTwoStrata(
  resA,
  resB,
  labelA = "A",
  labelB = "B",
  matrices = c("pcor", "copulaCor"),
  deltaThreshold = 0.05,
  edgeOnlyPcor = TRUE,
  edgeOnlyCor = FALSE
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

- matrices:

  Character vector of matrices to compare (`"pcor"`, `"copulaCor"`).

- deltaThreshold:

  Threshold for direction labels and bar charts.

- edgeOnlyPcor:

  If `TRUE`, pcor comparison uses adjacency edges only.

- edgeOnlyCor:

  If `FALSE` (default), copula cor compares all pairs.

## Value

List with comparison data frames keyed by matrix name.
