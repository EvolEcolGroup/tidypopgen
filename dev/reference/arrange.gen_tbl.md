# An arrange method for `gen_tibble` objects

An arrange method for `gen_tibble` objects

## Usage

``` r
# S3 method for class 'gen_tbl'
arrange(..., deparse.level = 1)
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
test_gt %>% arrange(id)
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 3
#>   id    population  genotypes
#>   <chr> <chr>      <vctr_SNP>
#> 1 a     pop1        [1,1,...]
#> 2 b     pop1        [2,1,...]
#> 3 c     pop2        [2,.,...]
#> 4 d     pop2        [1,0,...]
#> 5 e     pop1        [1,2,...]
#> 6 f     pop3        [0,0,...]
#> 7 g     pop3        [0,1,...]
```
