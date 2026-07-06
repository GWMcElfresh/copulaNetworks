# Plot partial correlation heatmap

Plot partial correlation heatmap

## Usage

``` r
PlotPcorHeatmap(
  result,
  title = "Partial Correlation Matrix (Glasso)",
  vars = NULL,
  zeroDiag = TRUE
)
```

## Arguments

- result:

  Output of
  [`FitStratumCopula()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitStratumCopula.md).

- title:

  Plot title.

- vars:

  Optional subset of variables to display.

- zeroDiag:

  If `TRUE`, zero the diagonal for display.

## Value

pheatmap grob.
