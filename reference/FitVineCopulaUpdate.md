# Fit a vine copula update on uniform pseudo-observations (Phase 2a)

Requires the optional **VineCopula** package.

## Usage

``` r
FitVineCopulaUpdate(
  uniformMatrix,
  factorGroups = NULL,
  selectionCrit = c("AIC", "BIC")
)
```

## Arguments

- uniformMatrix:

  Matrix of pseudo-observations in (0, 1) (n x d).

- factorGroups:

  Optional named integer vector of factor groups (metadata).

- selectionCrit:

  Selection criterion for pair-copula families.

## Value

List with `vineFit`, `logLik`, `n`, `factorGroups`, and `familySet`.
