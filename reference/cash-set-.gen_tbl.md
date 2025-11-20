# A \$ method for `gen_tibble` objects

A \$ method for `gen_tibble` objects

## Usage

``` r
# S3 method for class 'gen_tbl'
x$i <- value
```

## Arguments

- x:

  a gen_tibble

- i:

  column name

- value:

  a value to assign

## Value

a `gen_tibble`

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Add a new column
example_gt$region <- "East"

example_gt
#> # A gen_tibble: 6 loci
#> # A tibble:     7 × 4
#>   id    population  genotypes region
#>   <chr> <chr>      <vctr_SNP> <chr> 
#> 1 a     pop1        [1,1,...] East  
#> 2 b     pop1        [2,1,...] East  
#> 3 c     pop2        [2,.,...] East  
#> 4 d     pop2        [1,0,...] East  
#> 5 e     pop1        [1,2,...] East  
#> 6 f     pop3        [0,0,...] East  
#> 7 g     pop3        [0,1,...] East  
```
