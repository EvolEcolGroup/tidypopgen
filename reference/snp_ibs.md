# Compute the Identity by State Matrix for a bigSNP object

This function computes the IBS matrix.

## Usage

``` r
snp_ibs(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  type = c("proportion", "adjusted_counts", "raw_counts"),
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

- type:

  one of "proportion" (equivalent to "ibs" in PLINK), "adjusted_counts"
  ("distance" in PLINK), and "raw_counts" (the counts of identical
  alleles and non-missing alleles, from which the two other quantities
  are computed)

- block.size:

  maximum number of columns read at once. Note that, to optimise the
  speed of matrix operations, we have to store in memory 3 times the
  columns.

## Value

if as.counts = TRUE function returns a list of two
[bigstatsr::FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
matrices, one of counts of IBS by alleles (i.e. 2\*n loci), and one of
valid alleles (i.e. 2 \* n_loci - 2 \* missing_loci). If as.counts =
FALSE returns a single matrix of IBS proportions.

## Details

Note that monomorphic sites are currently counted. Should we filter them
beforehand? What does plink do?

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

X <- attr(example_gt$genotypes, "fbm")
snp_ibs(X)
#>           [,1] [,2]  [,3]      [,4]      [,5] [,6]  [,7]
#> [1,] 1.0000000 0.80 0.700 0.8333333 0.7500000 0.60 0.700
#> [2,] 0.8000000 1.00 0.750 0.8000000 0.5000000 0.50 0.500
#> [3,] 0.7000000 0.75 1.000 0.6000000 0.7000000 0.75 0.625
#> [4,] 0.8333333 0.80 0.600 1.0000000 0.5833333 0.70 0.500
#> [5,] 0.7500000 0.50 0.700 0.5833333 1.0000000 0.60 0.500
#> [6,] 0.6000000 0.50 0.750 0.7000000 0.6000000 1.00 0.750
#> [7,] 0.7000000 0.50 0.625 0.5000000 0.5000000 0.75 1.000

# Compute for individuals 1 to 5
snp_ibs(X, ind.row = 1:5, ind.col = 1:5)
#>      [,1]      [,2]      [,3]  [,4]  [,5]
#> [1,] 1.00 0.7500000 0.7500000 0.800 0.800
#> [2,] 0.75 1.0000000 0.8333333 0.750 0.500
#> [3,] 0.75 0.8333333 1.0000000 0.625 0.625
#> [4,] 0.80 0.7500000 0.6250000 1.000 0.600
#> [5,] 0.80 0.5000000 0.6250000 0.600 1.000

# Adjust block.size
snp_ibs(X, block.size = 2)
#>           [,1] [,2]  [,3]      [,4]      [,5] [,6]  [,7]
#> [1,] 1.0000000 0.80 0.700 0.8333333 0.7500000 0.60 0.700
#> [2,] 0.8000000 1.00 0.750 0.8000000 0.5000000 0.50 0.500
#> [3,] 0.7000000 0.75 1.000 0.6000000 0.7000000 0.75 0.625
#> [4,] 0.8333333 0.80 0.600 1.0000000 0.5833333 0.70 0.500
#> [5,] 0.7500000 0.50 0.700 0.5833333 1.0000000 0.60 0.500
#> [6,] 0.6000000 0.50 0.750 0.7000000 0.6000000 1.00 0.750
#> [7,] 0.7000000 0.50 0.625 0.5000000 0.5000000 0.75 1.000

# Change type
snp_ibs(X, type = "proportion")
#>           [,1] [,2]  [,3]      [,4]      [,5] [,6]  [,7]
#> [1,] 1.0000000 0.80 0.700 0.8333333 0.7500000 0.60 0.700
#> [2,] 0.8000000 1.00 0.750 0.8000000 0.5000000 0.50 0.500
#> [3,] 0.7000000 0.75 1.000 0.6000000 0.7000000 0.75 0.625
#> [4,] 0.8333333 0.80 0.600 1.0000000 0.5833333 0.70 0.500
#> [5,] 0.7500000 0.50 0.700 0.5833333 1.0000000 0.60 0.500
#> [6,] 0.6000000 0.50 0.750 0.7000000 0.6000000 1.00 0.750
#> [7,] 0.7000000 0.50 0.625 0.5000000 0.5000000 0.75 1.000
snp_ibs(X, type = "adjusted_counts")
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]  6.0  4.8 4.20  5.0  4.5  3.6 4.20
#> [2,]  4.8  6.0 4.50  4.8  3.0  3.0 3.00
#> [3,]  4.2  4.5 6.00  3.6  4.2  4.5 3.75
#> [4,]  5.0  4.8 3.60  6.0  3.5  4.2 3.00
#> [5,]  4.5  3.0 4.20  3.5  6.0  3.6 3.00
#> [6,]  3.6  3.0 4.50  4.2  3.6  6.0 4.50
#> [7,]  4.2  3.0 3.75  3.0  3.0  4.5 6.00
snp_ibs(X, type = "raw_counts")
#> $ibs
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 
#> $valid_n
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 
```
