# Tidyverse methods for gt objects

A filter method for `gen_tibble` objects

## Usage

``` r
# S3 method for class 'gen_tbl'
filter(..., deparse.level = 1)
```

## Arguments

- ...:

  a gen_tibble and a data.frame or tibble

- deparse.level:

  an integer controlling the construction of column names.

## Value

a `gen_tibble`

## Examples

``` r
test_gt <- load_example_gt("gen_tbl")
test_gt %>% filter(id %in% c("a", "c"))
#> # A gen_tibble: 6 loci
#> # A tibble:     2 × 3
#>   id    population  genotypes
#>   <chr> <chr>      <vctr_SNP>
#> 1 a     pop1        [1,1,...]
#> 2 c     pop2        [2,.,...]
```
