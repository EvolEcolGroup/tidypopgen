# Subsetting method for matrix_int_names

Subsetting method for matrix_int_names

Assignment method for matrix_int_names

## Usage

``` r
# S3 method for class 'matrix_int_names'
x[i, j, i_names = NULL, j_names = NULL, drop = TRUE]

# S3 method for class 'matrix_int_names'
x[i, j, i_names = NULL, j_names = NULL] <- value
```

## Arguments

- x:

  A matrix_int_names object

- i:

  Row indices (numeric), row names (character), or missing)

- j:

  Column indices (numeric), column names (character), or missing)

- i_names:

  Optional integer row names to subset by

- j_names:

  Optional integer column names to subset by

- drop:

  Logical indicating whether to drop dimensions

- value:

  A value to assign to the subsetted matrix

## Value

A subsetted matrix_int_names object

## See also

Other matrix_int_names_functions:
[`as.data.frame.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/as.data.frame.matrix_int_names.md),
[`col_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/col_names.md),
[`matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/matrix_int_names.md),
[`print.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/print.matrix_int_names.md),
[`row_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/row_names.md)

## Examples

``` r
my_mat <- matrix_int_names(matrix(1:6, nrow = 2),
  row_names = c(10L, 20L),
  col_names = c(100L, 200L, 300L)
)
# subset by position
my_mat[1, 2:3, drop = FALSE]
#> Matrix with integer rownames and integer colnames
#> Dimensions: 1 x 2 
#>    200 300
#> 10   3   5
# subset by integer names
my_mat[i_names = 20L, j_names = c(100L, 300L), drop = FALSE]
#> Matrix with integer rownames and integer colnames
#> Dimensions: 1 x 2 
#>    100 300
#> 20   2   6
# assign by integer names
my_mat[i_names = 10L, j_names = 200L] <- 0L
my_mat
#> Matrix with integer rownames and integer colnames
#> Dimensions: 2 x 3 
#>    100 200 300
#> 10   1   0   5
#> 20   2   4   6
```
