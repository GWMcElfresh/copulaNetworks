# Two-Phase Factor-Vine Copula Networks

## Motivation

Modern copula models use latent factors to handle high-dimensional
dependence. In a **two-phase** workflow we fit a parsimonious factor
copula on a **large** reference cohort (Phase 1), then update with a
graphical model or vine copula on a **small** new sample (Phase 2),
testing whether the update is consistent with the prior. See
`factor_vine_plan.md` in the package source for literature context.

## Setup

``` r

library(copulaNetworks)
```

Optional dependencies (install as needed):

| Package        | Role                               |
|----------------|------------------------------------|
| FactorCopula   | Phase 1 factor copula MLE          |
| VineCopula     | Phase 2 vine update                |
| cmdstanr       | Bayesian Gaussian update           |
| bridgesampling | Marginal likelihood / Bayes factor |

``` r

has_factor <- requireNamespace("FactorCopula", quietly = TRUE)
has_vine <- requireNamespace("VineCopula", quietly = TRUE)
has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
```

## Synthetic data

``` r

set.seed(42)
d <- 12
N <- 500
n <- 30
lambda <- runif(d, 0.35, 0.75)
sim_block <- function(n_obs) {
  f <- rnorm(n_obs)
  eps <- matrix(rnorm(n_obs * d), nrow = n_obs, ncol = d)
  factor_part <- matrix(lambda, nrow = n_obs, ncol = d, byrow = TRUE) * f
  x <- factor_part + sweep(eps, 2, sqrt(1 - lambda^2), `*`)
  colnames(x) <- paste0("node", seq_len(d))
  as.data.frame(x)
}
priorData <- sim_block(N)
updateData <- sim_block(n)
nodes <- paste0("node", 1:d)
dim(priorData)
#> [1] 500  12
dim(updateData)
#> [1] 30 12
```

## Marginal transforms

Phase 2 must use **prior-fitted** marginal transforms, not re-ranking on
small `n`.

``` r

marginalSpec <- FitMarginalSpec(priorData, nodeCols = nodes)
U_update <- ApplyMarginalSpec(updateData, marginalSpec)
range(U_update)
#> [1] 0.006986028 0.997005988
```

Contrast with re-ranking on the small sample (wrong for copula update):

``` r

U_wrong <- UniformMarginalTransform(as.matrix(updateData))
mean(abs(U_update - U_wrong))
#> [1] 0.04511515
```

Normal scores for graphical update:

``` r

Z_update <- qnorm(U_update)
head(Z_update[, 1:3])
#>           node1      node2      node3
#> [1,] -0.4170433 -1.0860766  0.8044144
#> [2,] -1.2999656 -0.5175240 -1.1229069
#> [3,] -2.2911267 -0.9053699 -2.2911267
#> [4,]  0.9129336 -0.3100972  0.5638528
#> [5,] -0.5405434  0.8044144 -0.4061522
#> [6,]  0.3953091  0.1002326 -1.1913074
```

## Phase 1: Factor structure and prior fit

``` r

prior_fit <- FitFactorCopulaPrior(
  priorData,
  nodeCols = nodes,
  nFactors = 1L,
  linkingCopula = "bvn"
)
CheckFactorStructure(prior_fit$impliedCor, nFactors = 1L)
prior_fit$loadings
```

Implied correlation heatmap:

``` r

if (has_factor) {
  PlotCopulaCorHeatmap(
    list(copulaCor = prior_fit$impliedCor, keptCols = nodes),
    title = "Phase 1 implied correlation"
  )
}
```

## Phase 2(a): Vine update

``` r

update_fit_both <- FitCopulaUpdate(
  updateData,
  priorFit = prior_fit,
  method = "both",
  nlambda = 20
)
if (!is.null(update_fit_both$vine)) {
  update_fit_both$vine$logLik
}
FactorGroupsFromLoadings(prior_fit$loadings, k = 2)
```

## Phase 2(b): Graphical update

``` r

update_fit <- FitCopulaUpdate(
  updateData,
  priorFit = prior_fit,
  method = "graphical",
  nlambda = 20
)
update_fit$graphical$lambdaOpt
```

Network plot:

``` r

if (has_factor && !is.null(update_fit$graphical)) {
  PlotCopulaNetwork(update_fit$graphical, title = "update")
}
```

## Prior vs update comparison

``` r

cmp <- ComparePriorToUpdate(prior_fit, update_fit)
head(cmp$pcor)
```

## Frequentist consistency test

``` r

gof <- TestPriorConsistency(
  prior_fit,
  updateData,
  nRep = 200,
  seed = 1
)
gof$pValue
hist(gof$nullDistribution, main = "Null distribution", xlab = "pairwise tau stat")
abline(v = gof$obsStat, col = "red", lwd = 2)
```

## Bayesian update

``` r

if (has_factor && has_cmdstanr) {
  bayes_fit <- FitBayesianFactorUpdate(
    prior_fit,
    updateData,
    chains = 2,
    iter = 500,
    computeBayesFactor = requireNamespace("bridgesampling", quietly = TRUE),
    seed = 1
  )
  head(bayes_fit$summary)
}
```

## End-to-end pipeline

``` r

res <- RunFactorVinePipeline(
  priorData = priorData,
  updateData = updateData,
  nodeCols = nodes,
  nFactors = 1L,
  phase2Method = "graphical",
  testMethod = "simulation",
  nRep = 100,
  nlambda = 20,
  seed = 1
)
res$consistencyTest$pValue
```

Checkpoints use
[`SaveCheckpoint()`](https://GWMcElfresh.github.io/copulaNetworks/reference/SaveCheckpoint.md)
/
[`LoadCheckpoint()`](https://GWMcElfresh.github.io/copulaNetworks/reference/LoadCheckpoint.md)
when `outDir` is set.

## Relation to stratified API

For **single-table stratified** workflows, use
[`RunCopulaPipeline()`](https://GWMcElfresh.github.io/copulaNetworks/reference/RunCopulaPipeline.md)
with
[`PrepareCopulaData()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PrepareCopulaData.md)
and stratum specs. Use
[`RunFactorVinePipeline()`](https://GWMcElfresh.github.io/copulaNetworks/reference/RunFactorVinePipeline.md)
when you have separate large prior and small update cohorts.

## Limitations

- Covariate-adjusted marginals: user supplies pre-adjusted data in v1
- Continuous margins only for
  [`FitFactorCopulaPrior()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitFactorCopulaPrior.md)
- Factor-informed vine truncation uses heuristic
  [`FactorGroupsFromLoadings()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FactorGroupsFromLoadings.md)
- Full factor-tree vines deferred

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] copulaNetworks_0.3.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3       cli_3.6.6         knitr_1.51        rlang_1.2.0      
#>  [5] xfun_0.59         otel_0.2.0        generics_0.1.4    textshaping_1.0.5
#>  [9] jsonlite_2.0.0    glue_1.8.1        htmltools_0.5.9   ragg_1.5.2       
#> [13] sass_0.4.10       rmarkdown_2.31    tibble_3.3.1      evaluate_1.0.5   
#> [17] jquerylib_0.1.4   fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5  
#> [21] compiler_4.6.1    dplyr_1.2.1       fs_2.1.0          pkgconfig_2.0.3  
#> [25] systemfonts_1.3.2 digest_0.6.39     R6_2.6.1          tidyselect_1.2.1 
#> [29] pillar_1.11.1     magrittr_2.0.5    bslib_0.11.0      tools_4.6.1      
#> [33] pkgdown_2.2.0     cachem_1.1.0      desc_1.4.3
```
