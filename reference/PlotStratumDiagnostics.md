# Plot single-stratum diagnostics (network + heatmaps)

Plot single-stratum diagnostics (network + heatmaps)

## Usage

``` r
PlotStratumDiagnostics(
  fitResult,
  stratumLabel = "stratum",
  outDir = NULL,
  minPcor = 0.01,
  nodeGroups = NULL,
  seed = 42,
  width = 10,
  height = 10,
  networkWidth = NULL,
  networkHeight = NULL,
  heatmapWidth = NULL,
  heatmapHeight = NULL,
  dpi = 150
)
```

## Arguments

- fitResult:

  Output of
  [`FitStratumCopula()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitStratumCopula.md).

- stratumLabel:

  Label used in plot titles and file names.

- outDir:

  Directory for saved PNG/PDF files. If `NULL`, plots are not saved.

- minPcor:

  Minimum \|partial correlation\| for network edges.

- nodeGroups:

  Optional node group mapping.

- seed:

  Layout seed.

- width:

  Default save width in inches (network and heatmaps).

- height:

  Default save height in inches (network and heatmaps).

- networkWidth:

  Optional network plot width override.

- networkHeight:

  Optional network plot height override.

- heatmapWidth:

  Optional heatmap width override.

- heatmapHeight:

  Optional heatmap height override.

- dpi:

  Resolution for saved PNG files.

## Value

List with ggplot/grob objects: `network`, `copulaCorHeatmap`,
`pcorHeatmap`.
