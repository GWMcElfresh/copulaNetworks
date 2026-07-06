# Fit a factor copula prior on a large reference sample (Phase 1)

Requires the optional **FactorCopula** package. Continuous margins only
in v1; supply covariate-adjusted residuals before calling.

## Usage

``` r
FitFactorCopulaPrior(
  data,
  nodeCols,
  nFactors = 1L,
  linkingCopula = "bvn",
  nQuad = 25L
)
```

## Arguments

- data:

  Data frame containing node columns.

- nodeCols:

  Character vector of variable names.

- nFactors:

  Number of latent factors (1 or 2).

- linkingCopula:

  Character vector of linking copula families per variable (e.g. `"bvn"`
  for Gaussian).

- nQuad:

  Number of quadrature points for latent integration.

## Value

List with `factorFit`, `loadings`, `logLik`, `nFactors`,
`linkingCopula`, `marginalSpec`, and `impliedCor`.
