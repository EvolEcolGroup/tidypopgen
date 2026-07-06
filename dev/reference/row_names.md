# Get or set row names for matrix_int_names

Returns integer names if present, otherwise character dimnames. Use
instead of [`rownames()`](https://rdrr.io/r/base/colnames.html) for
guaranteed dispatch.

## Usage

``` r
row_names(x)

row_names(x) <- value

# Default S3 method
row_names(x) <- value

# S3 method for class 'matrix_int_names'
row_names(x) <- value
```

## Arguments

- x:

  A matrix_int_names object

- value:

  Integer or character vector of row names, or NULL

## Value

Integer or character vector of row names, or NULL

The modified matrix_int_names object

## See also

Other matrix_int_names_functions:
[`[.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/subset.md),
[`as.data.frame.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/as.data.frame.matrix_int_names.md),
[`col_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/col_names.md),
[`matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/matrix_int_names.md),
[`print.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/print.matrix_int_names.md)
