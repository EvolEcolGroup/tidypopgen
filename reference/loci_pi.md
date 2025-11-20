# Estimate nucleotide diversity (pi) at each locus

Estimate nucleotide diversity (pi) at each locus, accounting for missing
values. This uses the formula: c_0 \* c_1 / (n \* (n-1) / 2)

## Usage

``` r
loci_pi(.x, .col = "genotypes", n_cores, block_size, type, ...)

# S3 method for class 'tbl_df'
loci_pi(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
loci_pi(
  .x,
  .col = "genotypes",
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

# S3 method for class 'grouped_df'
loci_pi(
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
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

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

a vector of frequencies, one per locus

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# For pi
example_gt %>% loci_pi()
#> # A tibble: 18 × 3
#>    loci  group value
#>    <chr> <chr> <dbl>
#>  1 rs1   pop1  0.533
#>  2 rs1   pop2  0.5  
#>  3 rs1   pop3  0    
#>  4 rs2   pop1  0.533
#>  5 rs2   pop2  0    
#>  6 rs2   pop3  0.5  
#>  7 rs3   pop1  0    
#>  8 rs3   pop2  0    
#>  9 rs3   pop3  0.5  
#> 10 rs4   pop1  0.667
#> 11 rs4   pop2  0.5  
#> 12 rs4   pop3  0    
#> 13 rs5   pop1  0.6  
#> 14 rs5   pop2  0.5  
#> 15 rs5   pop3  1    
#> 16 rs6   pop1  0.333
#> 17 rs6   pop2  0.5  
#> 18 rs6   pop3  1    

# For pi per locus per population
example_gt %>%
  group_by(population) %>%
  loci_pi()
#> # A tibble: 18 × 3
#>    loci  group value
#>    <chr> <chr> <dbl>
#>  1 rs1   pop1  0.533
#>  2 rs1   pop2  0.5  
#>  3 rs1   pop3  0    
#>  4 rs2   pop1  0.533
#>  5 rs2   pop2  0    
#>  6 rs2   pop3  0.5  
#>  7 rs3   pop1  0    
#>  8 rs3   pop2  0    
#>  9 rs3   pop3  0.5  
#> 10 rs4   pop1  0.667
#> 11 rs4   pop2  0.5  
#> 12 rs4   pop3  0    
#> 13 rs5   pop1  0.6  
#> 14 rs5   pop2  0.5  
#> 15 rs5   pop3  1    
#> 16 rs6   pop1  0.333
#> 17 rs6   pop2  0.5  
#> 18 rs6   pop3  1    
# alternatively, return a list of populations with their pi
example_gt %>%
  group_by(population) %>%
  loci_pi(type = "list")
#> [[1]]
#> [1] 0.5333333 0.5333333 0.0000000 0.6666667 0.6000000 0.3333333
#> 
#> [[2]]
#> [1] 0.5 0.0 0.0 0.5 0.5 0.5
#> 
#> [[3]]
#> [1] 0.0 0.5 0.5 0.0 1.0 1.0
#> 
# or a matrix with populations in columns and loci in rows
example_gt %>%
  group_by(population) %>%
  loci_pi(type = "matrix")
#>          pop1 pop2 pop3
#> rs1 0.5333333  0.5  0.0
#> rs2 0.5333333  0.0  0.5
#> rs3 0.0000000  0.0  0.5
#> rs4 0.6666667  0.5  0.0
#> rs5 0.6000000  0.5  1.0
#> rs6 0.3333333  0.5  1.0
# or within reframe (not recommended, as it much less efficient
# than using it directly as shown above)
example_gt %>%
  group_by(population) %>%
  reframe(pi = loci_pi(genotypes))
#> # A tibble: 18 × 2
#>    population    pi
#>    <chr>      <dbl>
#>  1 pop1       0.533
#>  2 pop1       0.533
#>  3 pop1       0    
#>  4 pop1       0.667
#>  5 pop1       0.6  
#>  6 pop1       0.333
#>  7 pop2       0.5  
#>  8 pop2       0    
#>  9 pop2       0    
#> 10 pop2       0.5  
#> 11 pop2       0.5  
#> 12 pop2       0.5  
#> 13 pop3       0    
#> 14 pop3       0.5  
#> 15 pop3       0.5  
#> 16 pop3       0    
#> 17 pop3       1    
#> 18 pop3       1    
```
