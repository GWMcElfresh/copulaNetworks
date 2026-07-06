# copulaNetworks

Stratified and two-phase copula network modeling for high-dimensional dependence.

## Features

- **Stratified workflow** (`RunCopulaPipeline`): nonparanormal graphical lasso networks per stratum
- **Two-phase factor-vine** (`RunFactorVinePipeline`): large-sample factor copula prior + small-sample graphical/vine update
- **Bayesian meta-analysis** (`RunMetaAnalysisPipeline`): multiple historical cohorts via power prior + small-sample update

## Documentation

- [Two-Phase Factor-Vine Vignette](https://GWMcElfresh.github.io/copulaNetworks/articles/two-phase-factor-vine.html)
- [Meta-Analysis Power Prior Vignette](https://GWMcElfresh.github.io/copulaNetworks/articles/meta-analysis-power-prior.html)

## Quick start

```r
library(copulaNetworks)
# See vignette and inst/examples/example_factor_vine.R
```
