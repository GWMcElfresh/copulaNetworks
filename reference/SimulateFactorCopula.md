# Simulate data from a fitted factor copula prior

Simulate data from a fitted factor copula prior

## Usage

``` r
SimulateFactorCopula(priorFit, nSim = NULL, nObs = NULL)
```

## Arguments

- priorFit:

  Output of
  [`FitFactorCopulaPrior()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitFactorCopulaPrior.md).

- nSim:

  Number of simulated observations.

- nObs:

  Alias for `nSim` (either may be used).

## Value

Numeric matrix (nSim x d) of simulated continuous scores on the
normal-copula scale.
