# Nonparanormal transformation (rank-based Gaussian copula marginals)

Maps each column to normal scores via ranks: \$\$\hat{Z}\_j =
\Phi^{-1}(\mathrm{rank}(X_j) / (n + 1))\$\$

## Usage

``` r
NonparanormalTransform(inputMatrix)
```

## Arguments

- inputMatrix:

  Numeric matrix (n x p).

## Value

Matrix of normal scores with column names preserved.
