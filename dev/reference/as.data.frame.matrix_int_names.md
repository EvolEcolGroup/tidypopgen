# Coercion method to data.frame for matrix_int_names

This method converts a matrix_int_names object to a data.frame,
preserving the integer row names. Otherwise, it works similarly to the
default as.data.frame.matrix method.

## Usage

``` r
# S3 method for class 'matrix_int_names'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A matrix_int_names object

- row.names:

  Either a character or integer (as both are valid for data.frame). If
  NULL, row names are inherited from the `matrix_int_names`

- optional:

  Ignored, included for compatibility with generic signature

- ...:

  Additional arguments passed to as.data.frame.matrix

## See also

Other matrix_int_names_functions:
[`[.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/subset.md),
[`col_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/col_names.md),
[`matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/matrix_int_names.md),
[`print.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/print.matrix_int_names.md),
[`row_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/row_names.md)

## Examples

``` r
# Create a matrix with integer row and column names
my_mat <- matrix_int_names(matrix(1:6, nrow = 2),
  row_names = c(10L, 20L),
  col_names = c(100L, 200L, 300L)
)
# Convert to data.frame
my_df <- as.data.frame(my_mat)
my_df
#>    100 200 300
#> 10   1   3   5
#> 20   2   4   6
```
