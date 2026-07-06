# Run the full stratified copula pipeline

Convenience wrapper that runs prepare, fit, plot, and optional
comparisons. Each step writes RDS checkpoints when `outDir` is provided.

## Usage

``` r
RunCopulaPipeline(
  data,
  idCols = character(0),
  strataCols = character(0),
  nodeCols = NULL,
  excludeCols = character(0),
  strataSpecs,
  comparePairs = list(),
  outDir = "checkpoints/copula_run",
  method = c("stars", "ebic"),
  nlambda = 40,
  starsThresh = 0.1,
  minN = 10,
  includeFull = FALSE,
  plotDiagnostics = TRUE,
  nodeGroups = NULL,
  width = 10,
  height = 10,
  dpi = 150,
  comparisonWidth = 16,
  comparisonHeight = 11,
  comparisonDpi = 150,
  ...
)
```

## Arguments

- data:

  Clean input data frame.

- idCols:

  Identifier columns.

- strataCols:

  Exogenous stratification columns.

- nodeCols:

  Node columns (NULL = auto-detect numeric).

- excludeCols:

  Columns to exclude from modeling.

- strataSpecs:

  Named list of stratum recipes for
  [`PrepareCopulaData()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PrepareCopulaData.md).

- comparePairs:

  Optional list of length-2 character vectors naming strata to compare.

- outDir:

  Output directory for checkpoints and figures.

- method:

  Lambda selection method.

- nlambda:

  Number of lambda values.

- starsThresh:

  StARS threshold.

- minN:

  Minimum observations per stratum.

- includeFull:

  Fit unstratified baseline.

- plotDiagnostics:

  If `TRUE`, save per-stratum diagnostic plots.

- nodeGroups:

  Optional node group mapping for plots.

- width:

  Save width in inches for stratum diagnostic plots.

- height:

  Save height in inches for stratum diagnostic plots.

- dpi:

  Save resolution for stratum diagnostic PNG files.

- comparisonWidth:

  Save width in inches for comparison plots.

- comparisonHeight:

  Save height in inches for comparison plots.

- comparisonDpi:

  Save resolution for comparison PNG files.

- ...:

  Additional arguments passed to
  [`PlotStratumDiagnostics()`](https://GWMcElfresh.github.io/copulaNetworks/reference/PlotStratumDiagnostics.md).

## Value

List with `preparedData`, `fitResults`, `plotArtifacts`, and
`comparisons`.
