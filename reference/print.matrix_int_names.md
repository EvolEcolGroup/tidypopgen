# Print method for matrix_int_names

Print method for matrix_int_names

## Usage

``` r
# S3 method for class 'matrix_int_names'
print(x, ..., n = 10L)
```

## Arguments

- x:

  A matrix_int_names object

- ...:

  Additional arguments passed to print

- n:

  Optional number of rows/columns to print (default 10)

## Value

Invisibly returns the original object

## See also

Other matrix_int_names_functions:
[`[.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/reference/subset.md),
[`as.data.frame.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/reference/as.data.frame.matrix_int_names.md),
[`col_names()`](https://evolecolgroup.github.io/tidypopgen/reference/col_names.md),
[`matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/reference/matrix_int_names.md),
[`row_names()`](https://evolecolgroup.github.io/tidypopgen/reference/row_names.md)

## Examples

``` r
my_mat <- matrix_int_names(matrix(1:6, nrow = 2),
  row_names = c(10L, 20L),
  col_names = c(100L, 200L, 300L)
)
print(my_mat)
#> Matrix with integer rownames and integer colnames
#> Dimensions: 2 x 3 
#>    100 200 300
#> 10   1   3   5
#> 20   2   4   6
print(my_mat, n = 2)
#> Matrix with integer rownames and integer colnames
#> Dimensions: 2 x 3 
#>    100 200
#> 10   1   3
#> 20   2   4
#> # 1 more columns. Use `print(n = ...)` to see more
```
