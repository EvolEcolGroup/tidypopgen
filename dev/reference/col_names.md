# Get or set column names for matrix_int_names

This function reliably returns integer names if present, otherwise
character dimnames. Use this instead of colnames() if you need
guaranteed dispatch.

## Usage

``` r
col_names(x)

# Default S3 method
col_names(x)

# S3 method for class 'matrix_int_names'
col_names(x)

col_names(x) <- value

# Default S3 method
col_names(x) <- value

# S3 method for class 'matrix_int_names'
col_names(x) <- value
```

## Arguments

- x:

  A matrix_int_names object

- value:

  a vector of integer or character column names to set, replacing the
  current one.

## Value

Integer or character vector of column names, or NULL

## See also

Other matrix_int_names_functions:
[`[.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/subset.md),
[`as.data.frame.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/as.data.frame.matrix_int_names.md),
[`matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/matrix_int_names.md),
[`print.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/print.matrix_int_names.md),
[`row_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/row_names.md)
