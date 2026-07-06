# Plot diagnostics for all fitted strata

Plot diagnostics for all fitted strata

## Usage

``` r
PlotAllStrata(
  fits,
  outDir = "figures",
  width = 10,
  height = 10,
  dpi = 150,
  ...
)
```

## Arguments

- fits:

  Output of
  [`FitCopulaStrata()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitCopulaStrata.md)
  or a named list of fit results.

- outDir:

  Base output directory. Each stratum gets a subfolder.

- width:

  Save width in inches (passed to
  [`PlotStratumDiagnostics()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PlotStratumDiagnostics.md)).

- height:

  Save height in inches (passed to
  [`PlotStratumDiagnostics()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PlotStratumDiagnostics.md)).

- dpi:

  Save resolution for PNG files (passed to
  [`PlotStratumDiagnostics()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PlotStratumDiagnostics.md)).

- ...:

  Additional arguments passed to
  [`PlotStratumDiagnostics()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PlotStratumDiagnostics.md).

## Value

Named list of per-stratum plot objects.
