# Compute the Identity by State Matrix for a `gen_tibble` object

This function computes the IBS matrix.

## Usage

``` r
pairwise_ibs(
  x,
  as_matrix = FALSE,
  type = c("proportion", "adjusted_counts", "raw_counts"),
  block_size = bigstatsr::block_size(nrow(x))
)
```

## Arguments

- x:

  a `gen_tibble` object.

- as_matrix:

  boolean, determining whether the results should be a square
  symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)

- type:

  one of "proportion" (equivalent to "ibs" in PLINK), "adjusted_counts"
  ("distance" in PLINK), and "raw_counts" (the counts of identical
  alleles and non-missing alleles, from which the two other quantities
  are computed)

- block_size:

  maximum number of loci read at once. More loci should improve speed,
  but will tax memory.

## Value

a
[bigstatsr::FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
of proportion or adjusted counts, or a list of two
[bigstatsr::FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
matrices, one of counts of IBS by alleles, and one of number of valid
alleles (i.e. 2*n_loci - 2*missing_loci)

## Details

Note that monomorphic sites are currently considered. Remove monomorphic
sites before running pairwise_king if this is a concern.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

pairwise_ibs(example_gt, type = "proportion")
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b     0.8  
#>  2 a     c     0.7  
#>  3 a     d     0.833
#>  4 a     e     0.75 
#>  5 a     f     0.6  
#>  6 a     g     0.7  
#>  7 b     c     0.75 
#>  8 b     d     0.8  
#>  9 b     e     0.5  
#> 10 b     f     0.5  
#> # ℹ 11 more rows

# Alternatively, return a matrix
pairwise_ibs(example_gt, type = "proportion", as_matrix = TRUE)
#>           a    b     c         d         e    f     g
#> a 1.0000000 0.80 0.700 0.8333333 0.7500000 0.60 0.700
#> b 0.8000000 1.00 0.750 0.8000000 0.5000000 0.50 0.500
#> c 0.7000000 0.75 1.000 0.6000000 0.7000000 0.75 0.625
#> d 0.8333333 0.80 0.600 1.0000000 0.5833333 0.70 0.500
#> e 0.7500000 0.50 0.700 0.5833333 1.0000000 0.60 0.500
#> f 0.6000000 0.50 0.750 0.7000000 0.6000000 1.00 0.750
#> g 0.7000000 0.50 0.625 0.5000000 0.5000000 0.75 1.000
#> attr(,"class")
#> [1] "pairwise_matrix" "matrix"          "array"          

# Adjust block_size
pairwise_ibs(example_gt, block_size = 2)
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b     0.8  
#>  2 a     c     0.7  
#>  3 a     d     0.833
#>  4 a     e     0.75 
#>  5 a     f     0.6  
#>  6 a     g     0.7  
#>  7 b     c     0.75 
#>  8 b     d     0.8  
#>  9 b     e     0.5  
#> 10 b     f     0.5  
#> # ℹ 11 more rows

# Change type
pairwise_ibs(example_gt, type = "adjusted_counts")
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b       4.8
#>  2 a     c       4.2
#>  3 a     d       5  
#>  4 a     e       4.5
#>  5 a     f       3.6
#>  6 a     g       4.2
#>  7 b     c       4.5
#>  8 b     d       4.8
#>  9 b     e       3  
#> 10 b     f       3  
#> # ℹ 11 more rows
pairwise_ibs(example_gt, type = "raw_counts")
#> $ibs
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 
#> $valid_n
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 
```
