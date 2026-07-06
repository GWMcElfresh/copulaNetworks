# Transform each column to uniform margins via ranks

Transform each column to uniform margins via ranks

## Usage

``` r
UniformMarginalTransform(inputMatrix, nRef = NULL)
```

## Arguments

- inputMatrix:

  Numeric matrix (n x p).

- nRef:

  Reference sample size for rank denominator (default nrow of matrix).

## Value

Matrix of pseudo-observations in (0, 1) with column names preserved.
