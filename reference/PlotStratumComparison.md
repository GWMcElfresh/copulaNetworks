# Plot stratum comparison diagnostics

Plot stratum comparison diagnostics

## Usage

``` r
PlotStratumComparison(
  cmp,
  outDir = NULL,
  labelThreshold = 0.05,
  maxLabels = 25,
  deltaThreshold = 0.05,
  width = 16,
  height = 11,
  scatterWidth = NULL,
  scatterHeight = NULL,
  barWidth = NULL,
  barHeight = NULL,
  dpi = 150
)
```

## Arguments

- cmp:

  Output of
  [`CompareTwoStrata()`](https://GWMcElfresh.github.io/copulaNetworks/reference/CompareTwoStrata.md).

- outDir:

  Optional directory to save PNG files.

- labelThreshold:

  Minimum \|delta\| for edge labels in scatter plots.

- maxLabels:

  Maximum number of extreme pairs to label in scatter plots.

- deltaThreshold:

  Minimum \|delta\| for bar charts.

- width:

  Default save width in inches (scatter and bar plots).

- height:

  Default save height in inches (scatter and bar plots).

- scatterWidth:

  Optional scatter plot width override.

- scatterHeight:

  Optional scatter plot height override.

- barWidth:

  Optional bar plot width override.

- barHeight:

  Optional bar plot height override.

- dpi:

  Resolution for saved PNG files.

## Value

List with `plots` keyed by matrix type (each containing `scatter` and
`bar`).
