# copulaNetworks

Stratified and two-phase copula network modeling for high-dimensional dependence.

## Features

- **Stratified workflow** (`RunCopulaPipeline`): nonparanormal graphical lasso networks per stratum
- **Two-phase factor-vine** (`RunFactorVinePipeline`): large-sample factor copula prior + small-sample graphical/vine update
- **Bayesian meta-analysis** (`RunMetaAnalysisPipeline`): multiple historical cohorts via power prior + small-sample update
- **Latent mixed models** (`FitCopulaLatentMixedModel`): random effects in the latent Gaussian layer via `latentFormula` (brms HMC or conjugate Gibbs); residual NIW / glasso for the copula

## Documentation

- [Stratified Nonparanormal Vignette](https://GWMcElfresh.github.io/copulaNetworks/articles/stratified-nonparanormal.html)
- [Two-Phase Factor-Vine Vignette](https://GWMcElfresh.github.io/copulaNetworks/articles/two-phase-factor-vine.html)
- [Meta-Analysis Power Prior Vignette](https://GWMcElfresh.github.io/copulaNetworks/articles/meta-analysis-power-prior.html)
- [Copula Mixed Models Vignette](https://GWMcElfresh.github.io/copulaNetworks/articles/copula-mixed-models.html)

## Quick start

```r
library(copulaNetworks)
# See vignettes and inst/examples/
```
