# Compute the Identity by State Matrix for a bigSNP object

This function computes the IBS matrix.

## Usage

``` r
snp_ibs(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  type = c("proportion", "adjusted_counts", "raw_counts"),
  block.size = bigstatsr::block_size(nrow(X)),
  n_cores = bigstatsr::nb_cores()
)
```

## Arguments

- X:

  a
  [bigstatsr::FBM.code256](https://privefl.github.io/bigstatsr/reference/FBM.code256-class.html)
  matrix (as found in the `genotypes` slot of a
  [bigsnpr::bigSNP](https://privefl.github.io/bigsnpr/reference/bigSNP-class.html)
  object). Genotypes are decoded with its `code256`, and values other
  than 0, 1 and 2 are treated as missing.

- ind.row:

  An optional vector of the row indices that are used. If not specified,
  all rows are used. Don't use negative indices.

- ind.col:

  An optional vector of the column indices that are used. If not
  specified, all columns are used. Don't use negative indices.

- type:

  one of:

  - "proportion" (the default): the proportion of alleles that are
    identical by state, a similarity (as `--distance ibs flat-missing`
    in PLINK, which writes it to a .mibs file).

  - "adjusted_counts": the allele count distance, a dissimilarity (as
    `--distance flat-missing` in PLINK, which writes it to a .dist
    file): the number of alleles that differ between two individuals,
    scaled up by the number of loci over the number of loci genotyped in
    both individuals, to correct for missing genotypes. It equals
    `2 * n_loci * (1 - proportion)`, where `n_loci` is the number of
    loci used, so it is 0 for identical individuals and at most
    `2 * n_loci`. Note that, without `flat-missing`, PLINK corrects for
    missing genotypes differently, which gives slightly different values
    when genotypes are missing.

  - "raw_counts": the counts of identical alleles and of non-missing
    alleles, from which the two other quantities are computed.

- block.size:

  maximum number of columns read at once. Genotypes are compared in a
  compact bit-packed form, so this has little effect on memory use.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

## Value

for "proportion" and "adjusted_counts", a matrix with one row and one
column per individual; for "raw_counts", a list of two
[bigstatsr::FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
matrices, `ibs` (the counts of alleles identical by state) and `valid_n`
(the counts of non-missing alleles, i.e. twice the number of loci
genotyped in both individuals).

## Details

Note that monomorphic sites are counted (as in PLINK). Remove
monomorphic sites beforehand if this is a concern.

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
#> [1,]  0.0  2.4  3.6  2.0  3.0  4.8  3.6
#> [2,]  2.4  0.0  3.0  2.4  6.0  6.0  6.0
#> [3,]  3.6  3.0  0.0  4.8  3.6  3.0  4.5
#> [4,]  2.0  2.4  4.8  0.0  5.0  3.6  6.0
#> [5,]  3.0  6.0  3.6  5.0  0.0  4.8  6.0
#> [6,]  4.8  6.0  3.0  3.6  4.8  0.0  3.0
#> [7,]  3.6  6.0  4.5  6.0  6.0  3.0  0.0
snp_ibs(X, type = "raw_counts")
#> $ibs
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 
#> $valid_n
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 
```
