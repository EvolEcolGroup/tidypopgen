# Matrix with integer row/column names class

Create a matrix with integer or character names for rows and columns.
Integer names are useful when dealing with very large numbers of rows or
columns, as the integer names can be stored more efficiently than
character names. Whilst this is possible for `data.frames`, the native
`matrix` class only allows for character row and column names via the
`dimnames` attribute. Use
[`row_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/row_names.md)
and
[`col_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/col_names.md)
to get or set the integer names, which are stored in special attributes
(`int_rownames` and `int_colnames`) on the matrix.

## Usage

``` r
matrix_int_names(data, row_names = NULL, col_names = NULL)
```

## Arguments

- data:

  Matrix data or object coercible to matrix

- row_names:

  Integer vector (stored as int_rownames) or character vector (stored as
  dimnames)

- col_names:

  Integer vector (stored as int_colnames) or character vector (stored as
  dimnames)

## Value

A `matrix_int_names` object

## Details

This class allows you to have integer names as attributes while still
being a matrix, and provides methods for getting and setting these
names. When you create a `matrix_int_names` object, you can provide
integer or character vectors for row and column names. If you provide
integer vectors, they will be stored in special attributes
(`int_rownames` and `int_colnames`) instead of the standard `dimnames`.
If you provide character vectors, they will be stored in the standard
`dimnames` as usual. You can also mix and match, having integer names
for rows and character names for columns, or vice versa. Note that,
since the row and column names are stored in special attributes, you
have to use
[`row_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/row_names.md)
and
[`col_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/col_names.md)
to get or set them, rather than
[`rownames()`](https://rdrr.io/r/base/colnames.html) and
[`colnames()`](https://rdrr.io/r/base/colnames.html), which will only
return character names if present.

## See also

Other matrix_int_names_functions:
[`[.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/subset.md),
[`as.data.frame.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/as.data.frame.matrix_int_names.md),
[`col_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/col_names.md),
[`print.matrix_int_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/print.matrix_int_names.md),
[`row_names()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/row_names.md)

## Examples

``` r
# Create a matrix with integer row and column names
my_mat <- matrix_int_names(matrix(1:6, nrow = 2),
  row_names = c(10L, 20L),
  col_names = c(100L, 200L, 300L)
)
row_names(my_mat) # returns integer row names
#> [1] 10 20
col_names(my_mat) # returns integer column names
#> [1] 100 200 300
```
