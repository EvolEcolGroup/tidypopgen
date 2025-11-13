# Estimate missingness at each locus

Estimate the rate of missingness at each locus. This function has an
efficient method to support grouped `gen_tibble` objects, which can
return a tidied tibble, a list, or a matrix.

## Usage

``` r
loci_missingness(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size,
  type,
  ...
)

# S3 method for class 'tbl_df'
loci_missingness(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
loci_missingness(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

# S3 method for class 'grouped_df'
loci_missingness(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotypes` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- as_counts:

  boolean defining whether the count of NAs (rather than the rate)
  should be returned. It defaults to FALSE (i.e. rates are returned by
  default).

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html)

- block_size:

  maximum number of loci read at once.

- type:

  type of object to return, if using grouped method. One of "tidy",
  "list", or "matrix". Default is "tidy".

- ...:

  other arguments passed to specific methods.

## Value

a vector of frequencies, one per locus

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# For missingness
example_gt %>% loci_missingness()
#> [1] 0.0000000 0.1428571 0.0000000 0.1428571 0.1428571 0.1428571

# For missingness per locus per population
example_gt %>%
  group_by(population) %>%
  loci_missingness()
#> # A tibble: 18 × 3
#>    loci  group value
#>    <chr> <chr> <dbl>
#>  1 rs1   pop1  0    
#>  2 rs1   pop2  0    
#>  3 rs1   pop3  0    
#>  4 rs2   pop1  0    
#>  5 rs2   pop2  0.5  
#>  6 rs2   pop3  0    
#>  7 rs3   pop1  0    
#>  8 rs3   pop2  0    
#>  9 rs3   pop3  0    
#> 10 rs4   pop1  0.333
#> 11 rs4   pop2  0    
#> 12 rs4   pop3  0    
#> 13 rs5   pop1  0    
#> 14 rs5   pop2  0    
#> 15 rs5   pop3  0.5  
#> 16 rs6   pop1  0    
#> 17 rs6   pop2  0    
#> 18 rs6   pop3  0.5  
# alternatively, return a list of populations with their missingness
example_gt %>%
  group_by(population) %>%
  loci_missingness(type = "list")
#> [[1]]
#> [1] 0.0000000 0.0000000 0.0000000 0.3333333 0.0000000 0.0000000
#> 
#> [[2]]
#> [1] 0.0 0.5 0.0 0.0 0.0 0.0
#> 
#> [[3]]
#> [1] 0.0 0.0 0.0 0.0 0.5 0.5
#> 
# or a matrix with populations in columns and loci in rows
example_gt %>%
  group_by(population) %>%
  loci_missingness(type = "matrix")
#>          pop1 pop2 pop3
#> rs1 0.0000000  0.0  0.0
#> rs2 0.0000000  0.5  0.0
#> rs3 0.0000000  0.0  0.0
#> rs4 0.3333333  0.0  0.0
#> rs5 0.0000000  0.0  0.5
#> rs6 0.0000000  0.0  0.5
# or within reframe (not recommended, as it much less efficient
# than using it directly as shown above)
example_gt %>%
  group_by(population) %>%
  reframe(missing = loci_missingness(genotypes))
#> # A tibble: 18 × 2
#>    population missing
#>    <chr>        <dbl>
#>  1 pop1         0    
#>  2 pop1         0    
#>  3 pop1         0    
#>  4 pop1         0.333
#>  5 pop1         0    
#>  6 pop1         0    
#>  7 pop2         0    
#>  8 pop2         0.5  
#>  9 pop2         0    
#> 10 pop2         0    
#> 11 pop2         0    
#> 12 pop2         0    
#> 13 pop3         0    
#> 14 pop3         0    
#> 15 pop3         0    
#> 16 pop3         0    
#> 17 pop3         0.5  
#> 18 pop3         0.5  
```
