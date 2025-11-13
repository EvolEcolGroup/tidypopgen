# Estimate allele frequencies at each locus

Allele frequencies can be estimates as minimum allele frequencies (MAF)
with `loci_maf()` or the frequency of the alternate allele (with
`loci_alt_freq()`). The latter are in line with the genotypes matrix
(e.g. as extracted by
[`show_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/show_loci.md)).
Most users will be in interested in the MAF, but the raw frequencies
might be useful when computing aggregated statistics. Both `loci_maf()`
and `loci_alt_freq()` have efficient methods to support grouped
`gen_tibble` objects. These can return a tidied tibble, a list, or a
matrix.

## Usage

``` r
loci_alt_freq(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores,
  block_size,
  type,
  ...
)

# S3 method for class 'tbl_df'
loci_alt_freq(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
loci_alt_freq(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

# S3 method for class 'grouped_df'
loci_alt_freq(
  .x,
  .col = "genotypes",
  as_counts = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  type = c("tidy", "list", "matrix"),
  ...
)

loci_maf(.x, .col = "genotypes", n_cores, block_size, type, ...)

# S3 method for class 'tbl_df'
loci_maf(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
loci_maf(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

# S3 method for class 'grouped_df'
loci_maf(
  .x,
  .col = "genotypes",
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

  boolean defining whether the count of alternate and valid (i.e. total
  number) alleles (rather than the frequencies) should be returned. It
  defaults to FALSE (i.e. frequencies are returned by default).

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html)

- block_size:

  maximum number of loci read at once.

- type:

  type of object to return, if using grouped method. One of "tidy",
  "list", or "matrix". Default is "tidy".

- ...:

  other arguments passed to specific methods, currently unused.

## Value

a vector of frequencies, one per locus, if `as_counts = FALSE`; else a
matrix of two columns, the count of alternate alleles and the count
valid alleles (i.e. the sum of alternate and reference)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# For alternate allele frequency
example_gt %>% loci_alt_freq()
#> [1] 0.50000000 0.41666667 0.07142857 0.25000000 0.41666667 0.25000000

# For alternate allele frequency per locus per population
example_gt %>%
  group_by(population) %>%
  loci_alt_freq()
#> # A tibble: 18 × 3
#>    loci  group value
#>    <chr> <chr> <dbl>
#>  1 rs1   pop1  0.667
#>  2 rs1   pop2  0.75 
#>  3 rs1   pop3  0    
#>  4 rs2   pop1  0.667
#>  5 rs2   pop2  0    
#>  6 rs2   pop3  0.25 
#>  7 rs3   pop1  0    
#>  8 rs3   pop2  0    
#>  9 rs3   pop3  0.25 
#> 10 rs4   pop1  0.5  
#> 11 rs4   pop2  0.25 
#> 12 rs4   pop3  0    
#> 13 rs5   pop1  0.5  
#> 14 rs5   pop2  0.25 
#> 15 rs5   pop3  0.5  
#> 16 rs6   pop1  0.167
#> 17 rs6   pop2  0.25 
#> 18 rs6   pop3  0.5  
# alternatively, return a list of populations with their frequencies
example_gt %>%
  group_by(population) %>%
  loci_alt_freq(type = "list")
#> [[1]]
#> [1] 0.6666667 0.6666667 0.0000000 0.5000000 0.5000000 0.1666667
#> 
#> [[2]]
#> [1] 0.75 0.00 0.00 0.25 0.25 0.25
#> 
#> [[3]]
#> [1] 0.00 0.25 0.25 0.00 0.50 0.50
#> 
# or a matrix with populations in columns and loci in rows
example_gt %>%
  group_by(population) %>%
  loci_alt_freq(type = "matrix")
#>          pop1 pop2 pop3
#> rs1 0.6666667 0.75 0.00
#> rs2 0.6666667 0.00 0.25
#> rs3 0.0000000 0.00 0.25
#> rs4 0.5000000 0.25 0.00
#> rs5 0.5000000 0.25 0.50
#> rs6 0.1666667 0.25 0.50
# or within reframe (not recommended, as it much less efficient
# than using it directly as shown above)
library(dplyr)
example_gt %>%
  group_by(population) %>%
  reframe(alt_freq = loci_alt_freq(genotypes))
#> # A tibble: 18 × 2
#>    population alt_freq
#>    <chr>         <dbl>
#>  1 pop1          0.667
#>  2 pop1          0.667
#>  3 pop1          0    
#>  4 pop1          0.5  
#>  5 pop1          0.5  
#>  6 pop1          0.167
#>  7 pop2          0.75 
#>  8 pop2          0    
#>  9 pop2          0    
#> 10 pop2          0.25 
#> 11 pop2          0.25 
#> 12 pop2          0.25 
#> 13 pop3          0    
#> 14 pop3          0.25 
#> 15 pop3          0.25 
#> 16 pop3          0    
#> 17 pop3          0.5  
#> 18 pop3          0.5  
# For MAF
example_gt %>% loci_maf()
#> [1] 0.50000000 0.41666667 0.07142857 0.25000000 0.41666667 0.25000000

# For minor allele frequency per locus per population
example_gt %>%
  group_by(population) %>%
  loci_maf()
#> # A tibble: 18 × 3
#>    loci  group value
#>    <chr> <chr> <dbl>
#>  1 rs1   pop1  0.333
#>  2 rs1   pop2  0.25 
#>  3 rs1   pop3  0    
#>  4 rs2   pop1  0.333
#>  5 rs2   pop2  0    
#>  6 rs2   pop3  0.25 
#>  7 rs3   pop1  0    
#>  8 rs3   pop2  0    
#>  9 rs3   pop3  0.25 
#> 10 rs4   pop1  0.5  
#> 11 rs4   pop2  0.25 
#> 12 rs4   pop3  0    
#> 13 rs5   pop1  0.5  
#> 14 rs5   pop2  0.25 
#> 15 rs5   pop3  0.5  
#> 16 rs6   pop1  0.167
#> 17 rs6   pop2  0.25 
#> 18 rs6   pop3  0.5  
# alternatively, return a list of populations with their frequencies
example_gt %>%
  group_by(population) %>%
  loci_maf(type = "list")
#> [[1]]
#> [1] 0.3333333 0.3333333 0.0000000 0.5000000 0.5000000 0.1666667
#> 
#> [[2]]
#> [1] 0.25 0.00 0.00 0.25 0.25 0.25
#> 
#> [[3]]
#> [1] 0.00 0.25 0.25 0.00 0.50 0.50
#> 
# or a matrix with populations in columns and loci in rows
example_gt %>%
  group_by(population) %>%
  loci_maf(type = "matrix")
#>          pop1 pop2 pop3
#> rs1 0.3333333 0.25 0.00
#> rs2 0.3333333 0.00 0.25
#> rs3 0.0000000 0.00 0.25
#> rs4 0.5000000 0.25 0.00
#> rs5 0.5000000 0.25 0.50
#> rs6 0.1666667 0.25 0.50
```
