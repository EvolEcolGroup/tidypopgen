# Compute the KING-robust Matrix for a bigSNP object

This function computes the KING-robust estimator of kinship,
reimplementing the KING algorithm of Manichaikul et al. (2010).

## Usage

``` r
snp_king(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  block.size = bigstatsr::block_size(nrow(X)) * 4
)
```

## Arguments

- X:

  a
  [bigstatsr::FBM.code256](https://privefl.github.io/bigstatsr/reference/FBM.code256-class.html)
  matrix (as found in the `genotypes` slot of a
  [bigsnpr::bigSNP](https://privefl.github.io/bigsnpr/reference/bigSNP-class.html)
  object).

- ind.row:

  An optional vector of the row indices that are used. If not specified,
  all rows are used. Don't use negative indices.

- ind.col:

  An optional vector of the column indices that are used. If not
  specified, all columns are used. Don't use negative indices.

- block.size:

  maximum number of columns read at once.

## Value

a square symmetrical matrix of relationship coefficients between
individuals

## References

Manichaikul, A. et al. (2010) Robust relationship inference in
genome-wide association studies. Bioinformatics, 26(22), 2867–2873.
https://doi.org/10.1093/bioinformatics/btq559.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

X <- attr(example_gt$genotypes, "fbm")
snp_king(X)
#>        [,1]  [,2]   [,3]   [,4]   [,5]  [,6]   [,7]
#> [1,]  0.500  0.00  0.125  0.250  0.250 -0.50  0.250
#> [2,]  0.000  0.50    NaN  0.000 -1.250 -1.00 -1.000
#> [3,]  0.125   NaN  0.500  0.000  0.125 -0.50 -0.750
#> [4,]  0.250  0.00  0.000  0.500 -0.625 -0.25 -0.125
#> [5,]  0.250 -1.25  0.125 -0.625  0.500 -1.00 -0.125
#> [6,] -0.500 -1.00 -0.500 -0.250 -1.000  0.50    NaN
#> [7,]  0.250 -1.00 -0.750 -0.125 -0.125   NaN  0.500

# Compute for individuals 1 to 5
snp_king(X, ind.row = 1:5, ind.col = 1:5)
#>      [,1] [,2]  [,3]  [,4]  [,5]
#> [1,] 0.50  0.0  0.00  0.25  0.25
#> [2,] 0.00  0.5   NaN  0.00 -1.00
#> [3,] 0.00  NaN  0.50 -0.25 -0.25
#> [4,] 0.25  0.0 -0.25  0.50 -0.50
#> [5,] 0.25 -1.0 -0.25 -0.50  0.50

# Adjust block size
snp_king(X, block.size = 2)
#>        [,1]  [,2]   [,3]   [,4]   [,5]  [,6]   [,7]
#> [1,]  0.500  0.00  0.125  0.250  0.250 -0.50  0.250
#> [2,]  0.000  0.50    NaN  0.000 -1.250 -1.00 -1.000
#> [3,]  0.125   NaN  0.500  0.000  0.125 -0.50 -0.750
#> [4,]  0.250  0.00  0.000  0.500 -0.625 -0.25 -0.125
#> [5,]  0.250 -1.25  0.125 -0.625  0.500 -1.00 -0.125
#> [6,] -0.500 -1.00 -0.500 -0.250 -1.000  0.50    NaN
#> [7,]  0.250 -1.00 -0.750 -0.125 -0.125   NaN  0.500
```
