# Extract edge list from a copula fit result

Extract edge list from a copula fit result

## Usage

``` r
ExtractEdges(result, label)
```

## Arguments

- result:

  Output of
  [`FitStratumCopula()`](https://GWMcElfresh.github.io/copulaNetworks/reference/FitStratumCopula.md).

- label:

  Stratum label.

## Value

Data frame with columns `from`, `to`, `pcor`, `stratum`, `edge`.
