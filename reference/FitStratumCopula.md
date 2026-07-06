# Fit a nonparanormal graphical model on one stratum

Fit a nonparanormal graphical model on one stratum

## Usage

``` r
FitStratumCopula(
  data,
  nodeCols,
  nlambda = 40,
  method = c("stars", "ebic"),
  starsThresh = 0.1,
  preTransformed = FALSE,
  quantileMatrix = NULL
)
```

## Arguments

- data:

  Data frame containing node columns.

- nodeCols:

  Character vector of endogenous variable names.

- nlambda:

  Number of lambda values for the glasso path.

- method:

  Lambda selection criterion: `"stars"` (StARS) or `"ebic"`.

- starsThresh:

  StARS stability threshold (used when `method = "stars"`).

- preTransformed:

  If `TRUE`, `data` columns are already normal scores (skip
  [`NonparanormalTransform()`](https://GWMcElfresh.github.io/copulaNetworks/reference/NonparanormalTransform.md)).
  Ignored when `quantileMatrix` is supplied.

- quantileMatrix:

  Optional pre-computed normal-score matrix (n x p).

## Value

List with correlation matrix, partial correlation matrix, selected
graph, etc. Returns `NULL` if fewer than 3 non-constant variables.
