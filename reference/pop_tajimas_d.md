# Estimate Tajima's D for the whole genome

Note that Tajima's D estimates from data that have been filtered or
ascertained can be difficult to interpret. This function should ideally
be used on sequence data prior to filtering.

## Usage

``` r
pop_tajimas_d(.x, n_cores, block_size, ...)

# S3 method for class 'tbl_df'
pop_tajimas_d(
  .x,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
pop_tajimas_d(
  .x,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

# S3 method for class 'grouped_df'
pop_tajimas_d(
  .x,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotypes` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html)

- block_size:

  maximum number of loci read at once.

- ...:

  other arguments passed to specific methods, currently unused.

## Value

A single numeric value (Tajima's D D) for the whole data set, `NA` when
the statistic is not defined. For grouped data a list of Tajima's D D
values (one per group) is returned.

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Compute Tajima's D
example_gt %>% pop_tajimas_d()
#> [[1]]
#> [1] 1.218829
#> 
#> [[2]]
#> [1] -0.780123
#> 
#> [[3]]
#> [1] 14.90782
#> 
```
