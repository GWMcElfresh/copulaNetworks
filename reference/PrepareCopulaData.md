# Prepare data and strata for copula modeling (Step 0)

Validates column roles, resolves node variables, and builds named strata
from declarative recipes. User supplies a clean data frame; no
imputation is performed.

## Usage

``` r
PrepareCopulaData(
  data,
  idCols = character(0),
  strataCols = character(0),
  nodeCols = NULL,
  excludeCols = character(0),
  strataSpecs = list(all = list()),
  outDir = NULL
)
```

## Arguments

- data:

  Clean input data frame.

- idCols:

  Character vector of identifier columns (excluded from model).

- strataCols:

  Character vector of exogenous stratification columns (excluded from
  model).

- nodeCols:

  Character vector of endogenous node columns. If `NULL`, all numeric
  columns not in `idCols`, `strataCols`, or `excludeCols` are used.

- excludeCols:

  Additional columns to exclude from the model.

- strataSpecs:

  Named list of stratum recipes passed to
  [`BuildStrata()`](https://GWMcElfresh.github.io/copulaNetworks/reference/BuildStrata.md).
  Each recipe may include `mutate`, `filter`, `group_by`, `stratumCol`,
  `nameSep`, `minN`. Keys become prefixes: stratum names are
  `"<spec_key>::<stratum_label>"`.

- outDir:

  Optional directory to save `prep.rds`.

## Value

List with elements `data`, `idCols`, `strataCols`, `nodeCols`, `strata`,
and `meta` (per-stratum row counts).
