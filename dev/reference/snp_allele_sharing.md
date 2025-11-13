# Compute the Pairwise Allele Sharing Matrix for a bigSNP object

This function computes the Allele Sharing matrix. Estimates Allele
Sharing (matching in `hierfstat`)) between pairs of individuals (for
each locus, gives 1 if the two individuals are homozygous for the same
allele, 0 if they are homozygous for a different allele, and 1/2 if at
least one individual is heterozygous. Matching is the average of these
0, 1/2 and 1s)

## Usage

``` r
snp_allele_sharing(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  block.size = bigstatsr::block_size(nrow(X))
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

  maximum number of columns read at once. Note that, to optimise the
  speed of matrix operations, we have to store in memory 3 times the
  columns.

## Value

a matrix of allele sharing between all pairs of individuals

## See also

[`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_allele_sharing.md)
[`hierfstat::matching()`](https://rdrr.io/pkg/hierfstat/man/matching.html)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

X <- attr(example_gt$genotypes, "fbm")
snp_allele_sharing(X)
#>           [,1]  [,2]  [,3]      [,4]      [,5]  [,6]  [,7]
#> [1,] 0.6666667 0.700 0.600 0.6666667 0.5833333 0.600 0.500
#> [2,] 0.7000000 0.900 0.750 0.8000000 0.5000000 0.500 0.375
#> [3,] 0.6000000 0.750 0.800 0.6000000 0.6000000 0.625 0.500
#> [4,] 0.6666667 0.800 0.600 0.8333333 0.4166667 0.700 0.500
#> [5,] 0.5833333 0.500 0.600 0.4166667 0.7500000 0.500 0.500
#> [6,] 0.6000000 0.500 0.625 0.7000000 0.5000000 0.900 0.750
#> [7,] 0.5000000 0.375 0.500 0.5000000 0.5000000 0.750 0.700

# Compute for individuals 1 to 5
snp_allele_sharing(X, ind.row = 1:5, ind.col = 1:5)
#>       [,1]      [,2]      [,3]  [,4]  [,5]
#> [1,] 0.600 0.6250000 0.6250000 0.600 0.600
#> [2,] 0.625 0.8750000 0.8333333 0.750 0.500
#> [3,] 0.625 0.8333333 0.8750000 0.625 0.625
#> [4,] 0.600 0.7500000 0.6250000 0.800 0.400
#> [5,] 0.600 0.5000000 0.6250000 0.400 0.800

# Adjust block size
snp_allele_sharing(X, block.size = 2)
#>           [,1]  [,2]  [,3]      [,4]      [,5]  [,6]  [,7]
#> [1,] 0.6666667 0.700 0.600 0.6666667 0.5833333 0.600 0.500
#> [2,] 0.7000000 0.900 0.750 0.8000000 0.5000000 0.500 0.375
#> [3,] 0.6000000 0.750 0.800 0.6000000 0.6000000 0.625 0.500
#> [4,] 0.6666667 0.800 0.600 0.8333333 0.4166667 0.700 0.500
#> [5,] 0.5833333 0.500 0.600 0.4166667 0.7500000 0.500 0.500
#> [6,] 0.6000000 0.500 0.625 0.7000000 0.5000000 0.900 0.750
#> [7,] 0.5000000 0.375 0.500 0.5000000 0.5000000 0.750 0.700
```
