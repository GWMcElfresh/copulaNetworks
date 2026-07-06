# Resolve node columns from data and role specifications

Resolve node columns from data and role specifications

## Usage

``` r
resolve_node_cols(data, nodeCols, idCols, strataCols, excludeCols)
```

## Arguments

- data:

  Input data frame.

- nodeCols:

  Explicit node columns, or NULL to auto-detect numeric columns.

- idCols:

  ID columns to exclude.

- strataCols:

  Stratification columns to exclude.

- excludeCols:

  Additional columns to exclude.

## Value

Character vector of node column names.
