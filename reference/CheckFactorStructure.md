# Check whether a correlation matrix is consistent with a low-rank factor structure

Check whether a correlation matrix is consistent with a low-rank factor
structure

## Usage

``` r
CheckFactorStructure(corMatrix, nFactors = 1L)
```

## Arguments

- corMatrix:

  Symmetric correlation matrix (d x d).

- nFactors:

  Number of factors to assess (default 1).

## Value

List with eigenvalues, proportion of variance explained, and a warning
flag when the first `nFactors` explain less than 50%.
