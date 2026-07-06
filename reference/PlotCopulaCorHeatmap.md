# Plot copula correlation heatmap

Plot copula correlation heatmap

## Usage

``` r
PlotCopulaCorHeatmap(result, title = "Copula Correlation Matrix", vars = NULL)
```

## Arguments

- result:

  Output of
  [`FitStratumCopula()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitStratumCopula.md).

- title:

  Plot title.

- vars:

  Optional subset of variables to display.

## Value

pheatmap grob.
