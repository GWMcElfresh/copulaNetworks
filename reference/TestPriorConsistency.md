# Test whether update data are consistent with the Phase 1 factor prior

Simulation-based goodness-of-fit: compares a pairwise Kendall tau
statistic on the real update sample to a null distribution from prior
replicates.

## Usage

``` r
TestPriorConsistency(
  priorFit,
  updateData,
  nRep = 500L,
  statistic = c("pairwiseTau"),
  seed = NULL
)
```

## Arguments

- priorFit:

  Output of
  [`FitFactorCopulaPrior()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitFactorCopulaPrior.md).

- updateData:

  Data frame with node columns for the update cohort.

- nRep:

  Number of simulation replicates (default 500).

- statistic:

  Test statistic: `"pairwiseTau"` (default).

- seed:

  Optional random seed for reproducibility.

## Value

List with `obsStat`, `pValue`, `nullDistribution`, and `nRep`.
