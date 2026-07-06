# Fit copula models across prepared strata

Fit copula models across prepared strata

## Usage

``` r
FitCopulaStrata(
  prep,
  method = c("stars", "ebic"),
  nlambda = 40,
  starsThresh = 0.1,
  minN = 10,
  includeFull = FALSE,
  outDir = NULL
)
```

## Arguments

- prep:

  Output of
  [`PrepareCopulaData()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PrepareCopulaData.md).

- method:

  Lambda selection criterion.

- nlambda:

  Number of lambda values.

- starsThresh:

  StARS threshold.

- minN:

  Minimum observations required per stratum.

- includeFull:

  If `TRUE`, also fit on the full (unstratified) dataset.

- outDir:

  Optional directory to save `fits.rds`.

## Value

Named list of fit results keyed by stratum name.
