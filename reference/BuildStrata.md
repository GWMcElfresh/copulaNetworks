# Build named strata from a declarative recipe

Build named strata from a declarative recipe

## Usage

``` r
BuildStrata(data, spec, specName = "strata")
```

## Arguments

- data:

  Input data frame.

- spec:

  Named list with optional elements:

  `mutate`

  :   Quoted expression evaluated with
      [`rlang::eval_tidy()`](https://rlang.r-lib.org/reference/eval_tidy.html).

  `filter`

  :   Quoted expression for row filtering.

  `group_by`

  :   Character vector of grouping columns, or a single column name.

  `stratumCol`

  :   Pre-built stratum column (skips `group_by` split).

  `nameSep`

  :   Separator when joining multiple `group_by` columns (default
      `" | "`).

  `minN`

  :   Minimum rows per stratum (default 1).

- specName:

  Name prefix for strata (used in messages).

## Value

Named list of stratum data frames.
