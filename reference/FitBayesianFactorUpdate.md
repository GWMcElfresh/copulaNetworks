# Bayesian Gaussian copula update with priors centered at Phase 1

Requires optional **cmdstanr** (and **bridgesampling** for Bayes
factors). ponytail: single-factor Gaussian copula via multivariate
normal on scores.

## Usage

``` r
FitBayesianFactorUpdate(
  priorFit,
  updateData,
  chains = 2L,
  iter = 1000L,
  computeBayesFactor = FALSE,
  seed = NULL
)
```

## Arguments

- priorFit:

  Output of
  [`FitFactorCopulaPrior()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitFactorCopulaPrior.md).

- updateData:

  Data frame with node columns for the update cohort.

- chains:

  Number of MCMC chains.

- iter:

  Total iterations per chain (warmup + sampling; split evenly).

- computeBayesFactor:

  If `TRUE` and bridgesampling is available, compute a
  marginal-likelihood estimate via bridge sampling on the fitted model.

- seed:

  Optional random seed.

## Value

List with `fit` (`CmdStanMCMC`), `summary`, and optional `bayesFactor`.
