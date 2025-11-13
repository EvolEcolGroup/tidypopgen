# Compute the Pairwise Allele Sharing Matrix for a `gen_tibble` object

This function computes the Allele Sharing matrix. Estimates Allele
Sharing (equivalent to the quantity estimated by
[`hierfstat::matching()`](https://rdrr.io/pkg/hierfstat/man/matching.html))
between pairs of individuals (for each locus, gives 1 if the two
individuals are homozygous for the same allele, 0 if they are homozygous
for a different allele, and 1/2 if at least one individual is
heterozygous. Matching is the average of these 0, 1/2 and 1s)

## Usage

``` r
pairwise_allele_sharing(
  x,
  as_matrix = FALSE,
  block_size = bigstatsr::block_size(nrow(x))
)
```

## Arguments

- x:

  a `gen_tibble` object.

- as_matrix:

  boolean, determining whether the results should be a square
  symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)

- block_size:

  maximum number of loci read at once. More loci should improve speed,
  but will tax memory.

## Value

a matrix of allele sharing between all pairs of individuals

## See also

[`hierfstat::matching()`](https://rdrr.io/pkg/hierfstat/man/matching.html)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Compute allele sharing between individuals
example_gt %>% pairwise_allele_sharing(as_matrix = FALSE)
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b     0.7  
#>  2 a     c     0.6  
#>  3 a     d     0.667
#>  4 a     e     0.583
#>  5 a     f     0.6  
#>  6 a     g     0.5  
#>  7 b     c     0.75 
#>  8 b     d     0.8  
#>  9 b     e     0.5  
#> 10 b     f     0.5  
#> # ℹ 11 more rows

# Alternatively, return as a tibble
example_gt %>% pairwise_allele_sharing(as_matrix = TRUE)
#>           a     b     c         d         e     f     g
#> a 0.6666667 0.700 0.600 0.6666667 0.5833333 0.600 0.500
#> b 0.7000000 0.900 0.750 0.8000000 0.5000000 0.500 0.375
#> c 0.6000000 0.750 0.800 0.6000000 0.6000000 0.625 0.500
#> d 0.6666667 0.800 0.600 0.8333333 0.4166667 0.700 0.500
#> e 0.5833333 0.500 0.600 0.4166667 0.7500000 0.500 0.500
#> f 0.6000000 0.500 0.625 0.7000000 0.5000000 0.900 0.750
#> g 0.5000000 0.375 0.500 0.5000000 0.5000000 0.750 0.700
```
