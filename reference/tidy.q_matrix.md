# Tidy a Q matrix

Takes a `q_matrix` object, which is a matrix, and returns a tidied
tibble.

## Usage

``` r
# S3 method for class 'q_matrix'
tidy(x, data, ...)
```

## Arguments

- x:

  A Q matrix object (as returned by
  [`q_matrix`](https://evolecolgroup.github.io/tidypopgen/reference/q_matrix.md)).

- data:

  An associated tibble (e.g. a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)),
  with the individuals in the same order as the data used to generate
  the Q matrix

- ...:

  not currently used

## Value

A tidied tibble containing columns:

- `row`:

  ID of the original observation (i.e. rowname from original data).

- `Q`:

  Integer indicating a Q component.

- `value`:

  The proportion for that particular Q value.

## Examples

``` r
# run the example only if we have the package installed
if (requireNamespace("LEA", quietly = TRUE)) {
  example_gt <- load_example_gt("gen_tbl")

  # Create a gt_admix object
  admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")

  # Extract a Q matrix
  q_mat_k3 <- get_q_matrix(admix_obj, k = 3, run = 1)

  tidy(q_mat_k3, data = example_gt)
}
#> # A tibble: 21 × 3
#>    id    q     percentage
#>    <chr> <chr>      <dbl>
#>  1 a     .Q1    1.000    
#>  2 a     .Q2    0.0001000
#>  3 a     .Q3    0.0001000
#>  4 b     .Q1    1.000    
#>  5 b     .Q2    0.0001000
#>  6 b     .Q3    0.0001000
#>  7 c     .Q1    0.0001000
#>  8 c     .Q2    1.000    
#>  9 c     .Q3    0.0001000
#> 10 d     .Q1    1.000    
#> # ℹ 11 more rows
```
