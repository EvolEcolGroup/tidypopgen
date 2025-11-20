# Compute the KING-robust Matrix for a `gen_tibble` object

This function computes the KING-robust estimator of kinship,
reimplementing the KING algorithm of Manichaikul et al. (2010).

## Usage

``` r
pairwise_king(
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

a square symmetrical matrix of relationship coefficients between
individuals if `as_matrix` is TRUE, or a tidied tibble of coefficients
if `as_matrix` is FALSE.

## References

Manichaikul, A. et al. (2010) Robust relationship inference in
genome-wide association studies. Bioinformatics, 26(22), 2867–2873.
https://doi.org/10.1093/bioinformatics/btq559.

Note that monomorphic sites are currently considered. Remove monomorphic
sites before running pairwise_king if this is a concern.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Compute the KING-robust matrix
pairwise_king(example_gt, as_matrix = TRUE)
#>        a     b      c      d      e     f      g
#> a  0.500  0.00  0.125  0.250  0.250 -0.50  0.250
#> b  0.000  0.50    NaN  0.000 -1.250 -1.00 -1.000
#> c  0.125   NaN  0.500  0.000  0.125 -0.50 -0.750
#> d  0.250  0.00  0.000  0.500 -0.625 -0.25 -0.125
#> e  0.250 -1.25  0.125 -0.625  0.500 -1.00 -0.125
#> f -0.500 -1.00 -0.500 -0.250 -1.000  0.50    NaN
#> g  0.250 -1.00 -0.750 -0.125 -0.125   NaN  0.500

# Or return a tidy tibble
pairwise_king(example_gt, as_matrix = FALSE)
#> # A tibble: 21 × 3
#>    item1 item2   value
#>    <chr> <chr>   <dbl>
#>  1 a     b       0    
#>  2 a     c       0.125
#>  3 a     d       0.25 
#>  4 a     e       0.25 
#>  5 a     f      -0.5  
#>  6 a     g       0.25 
#>  7 b     c     NaN    
#>  8 b     d       0    
#>  9 b     e      -1.25 
#> 10 b     f      -1    
#> # ℹ 11 more rows

# Adjust block_size
pairwise_king(example_gt, block_size = 2)
#> # A tibble: 21 × 3
#>    item1 item2   value
#>    <chr> <chr>   <dbl>
#>  1 a     b       0    
#>  2 a     c       0.125
#>  3 a     d       0.25 
#>  4 a     e       0.25 
#>  5 a     f      -0.5  
#>  6 a     g       0.25 
#>  7 b     c     NaN    
#>  8 b     d       0    
#>  9 b     e      -1.25 
#> 10 b     f      -1    
#> # ℹ 11 more rows
```
