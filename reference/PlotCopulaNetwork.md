# Plot a copula network using ggraph

Edge width and alpha encode \|partial correlation\|; edge colour encodes
sign.

## Usage

``` r
PlotCopulaNetwork(
  result,
  title = "",
  seed = 42,
  minPcor = 0.01,
  nodeGroups = NULL,
  printPlot = TRUE
)
```

## Arguments

- result:

  Output of
  [`FitStratumCopula()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitStratumCopula.md).

- title:

  Plot title.

- seed:

  Random seed for FR layout.

- minPcor:

  Minimum \|partial correlation\| to display an edge.

- nodeGroups:

  Optional group mapping (named vector or function).

- printPlot:

  If `TRUE`, print the plot before returning.

## Value

ggplot object, or `NULL` if nothing to plot.
