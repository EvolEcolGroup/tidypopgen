# Filter individuals based on a relationship threshold

This function takes a matrix of x by y individuals containing
relatedness coefficients and returns the maximum set of individuals that
contains no relationships above the given threshold.

## Usage

``` r
filter_high_relatedness(
  matrix,
  .x = NULL,
  kings_threshold = NULL,
  verbose = FALSE
)
```

## Arguments

- matrix:

  a square symmetric matrix of individuals containing relationship
  coefficients

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object

- kings_threshold:

  a threshold over which

- verbose:

  boolean whether to report to screen

## Value

a list where '1' is individual ID's to retain, '2' is individual ID's to
remove, and '3' is a boolean where individuals to keep are TRUE and
individuals to remove are FALSE

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Calculate relationship matrix
king_matrix <- example_gt %>% pairwise_king(as_matrix = TRUE)

# Filter individuals with threshold above 0.2
filter_high_relatedness(king_matrix, example_gt, kings_threshold = 0.2)
#> [[1]]
#> [1] "c" "d"
#> 
#> [[2]]
#> [1] "a" "b" "e" "f" "g"
#> 
#> [[3]]
#> [1] FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE
#> 
```
