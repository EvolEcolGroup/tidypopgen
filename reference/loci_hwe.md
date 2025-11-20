# Test Hardy-Weinberg equilibrium at each locus

Return the p-value from an exact test of HWE.

## Usage

``` r
loci_hwe(.x, .col = "genotypes", ...)

# S3 method for class 'tbl_df'
loci_hwe(.x, .col = "genotypes", mid_p = TRUE, ...)

# S3 method for class 'vctrs_bigSNP'
loci_hwe(.x, .col = "genotypes", mid_p = TRUE, ...)

# S3 method for class 'grouped_df'
loci_hwe(
  .x,
  .col = "genotypes",
  mid_p = TRUE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotypes` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- ...:

  not used.

- mid_p:

  boolean on whether the mid-p value should be computed. Default is
  TRUE, as in PLINK.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html)

- block_size:

  maximum number of loci read at once.

- type:

  type of object to return, if using grouped method. One of "tidy",
  "list", or "matrix". Default is "tidy".

## Value

a vector of probabilities from HWE exact test, one per locus

## Details

This function uses the original C++ algorithm from PLINK 1.90.

## Author

the C++ algorithm was written by Christopher Chang for PLINK 1.90, based
on original code by Jan Wigginton (the code was released under GPL3).

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# For HWE
example_gt %>% loci_hwe()
#> [1] 0.7552448 0.6969697 0.5000000 0.6363636 0.6969697 0.6363636

# For loci_hwe per locus per population, use reframe
example_gt %>%
  group_by(population) %>%
  reframe(loci_hwe = loci_hwe(genotypes))
#> # A tibble: 18 × 2
#>    population loci_hwe
#>    <chr>         <dbl>
#>  1 pop1          0.6  
#>  2 pop1          0.6  
#>  3 pop1          0.5  
#>  4 pop1          0.667
#>  5 pop1          0.7  
#>  6 pop1          0.5  
#>  7 pop2          0.5  
#>  8 pop2          0.5  
#>  9 pop2          0.5  
#> 10 pop2          0.5  
#> 11 pop2          0.5  
#> 12 pop2          0.5  
#> 13 pop3          0.5  
#> 14 pop3          0.5  
#> 15 pop3          0.5  
#> 16 pop3          0.5  
#> 17 pop3          0.5  
#> 18 pop3          0.5  
```
