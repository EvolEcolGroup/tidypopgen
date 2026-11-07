# Compute the joint site frequency spectrum of each pair of populations

Compute the two-dimensional (joint) site frequency spectrum (SFS) of
each pair of populations, i.e. the number of loci with \\i\\ copies of
the derived (or minor) allele in the first population and \\j\\ in the
second. Missing genotypes are dealt with as in `easySFS`, by projecting
each population down to `n_proj` allele copies (see
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)
for details). Each spectrum only uses the loci that have at least
`n_proj` called allele copies in both its populations, so its marginal
spectra are not, in general, equal to the spectra computed by
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md).
A folded spectrum is folded on the total number of copies across the two
populations (as in `dadi`), so that it counts copies of the allele that
is the minor one across the pair.

## Usage

``` r
pairwise_pop_sfs(.x, ...)

# S3 method for class 'gen_tbl'
pairwise_pop_sfs(
  .x,
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)
```

## Arguments

- .x:

  a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  (as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)),
  with at least two groups.

- ...:

  other arguments passed to specific methods, currently unused.

- n_proj:

  the number of allele copies to project each population to: either a
  single value used for all populations, a vector with one value per
  population (in the order of the groups), or a vector named with the
  group labels. If `NULL` (the default), each population is projected to
  the number of allele copies that maximises its number of segregating
  sites (see
  [`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md)).

- folded:

  boolean, whether to fold the spectrum (the default, as the ancestral
  allele is not usually known). If `FALSE`, the spectrum is polarised
  with `ancestral`.

- ancestral:

  the ancestral allele of each locus, used to polarise unfolded spectra:
  either the name of a column of the loci table (see
  [`show_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.md)),
  or a character vector with one allele per locus (in the order of the
  loci table). Loci whose ancestral allele is missing (`NA`) or is
  neither of their alleles are dropped. If `NULL` (the default), the
  reference allele is taken to be ancestral, as in `easySFS`. Ignored,
  with a warning, for folded spectra.

- type:

  the type of object to return, either "tidy" (the default) or "array".

- fastest:

  the order of the cells of each spectrum in the tidy output: "last"
  (the default) keeps the order of the populations, so that the rows are
  sorted by the number of allele copies in the first population, then in
  the second (the number of copies in the second population changes
  fastest, as in the files written by `dadi` and `easySFS`); "first"
  makes the number of copies in the first population change fastest.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- block_size:

  maximum number of loci read at once.

## Value

If `type = "tidy"`, a tibble with one row per cell of the joint spectrum
of each pair of populations, with the columns:

- the groups of the pair (e.g. `population_1` and `population_2`; one
  column per grouping variable for each population of the pair, if the
  data are grouped by more than one variable);

- `allele_count_1` and `allele_count_2`: the number of copies of the
  derived (or minor) allele in the first and second population (see
  `allele_count` in
  [`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md));

- `n_sites`: the (projected) number of sites;

- `masked`: whether the cell is masked (see
  [`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)).

The projection of each population is stored in the attribute `n_proj`.
The tibble has class `pairwise_sfs_tbl`, and can be plotted with
[autoplot()](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_sfs.md).
If `type = "array"`, a list of `sfs` objects (one per pair, named with
the labels of the two groups pasted with "."), each an array with the
first population in the rows and the second in the columns (see
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)).
Pairs are in the order given by
[`utils::combn()`](https://rdrr.io/r/utils/combn.html).

## See also

[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
[`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/nwise_pop_sfs.md),
[`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

example_gt %>% pairwise_pop_sfs(n_proj = 2)
#> # A tibble: 27 × 6
#>    population_1 population_2 allele_count_1 allele_count_2 n_sites masked
#>    <chr>        <chr>                 <int>          <int>   <dbl> <lgl> 
#>  1 pop1         pop2                      0              0   1.78  TRUE  
#>  2 pop1         pop2                      0              1   0.933 FALSE 
#>  3 pop1         pop2                      0              2   0.308 FALSE 
#>  4 pop1         pop2                      1              0   1.6   FALSE 
#>  5 pop1         pop2                      1              1   1.07  FALSE 
#>  6 pop1         pop2                      1              2   0     TRUE  
#>  7 pop1         pop2                      2              0   0.308 FALSE 
#>  8 pop1         pop2                      2              1   0     TRUE  
#>  9 pop1         pop2                      2              2   0     TRUE  
#> 10 pop1         pop3                      0              0   0.767 TRUE  
#> # ℹ 17 more rows

# as arrays
example_gt %>% pairwise_pop_sfs(n_proj = 2, type = "array")
#> $pop1.pop2
#> # A folded site frequency spectrum of 2 populations
#> # Allele copies: pop1 = 2, pop2 = 2
#> # Loci used: 6 (dropped for missing data: 0)
#>     pop2
#> pop1         0         1         2
#>    0 1.7833333 0.9333333 0.3083333
#>    1 1.6000000 1.0666667 0.0000000
#>    2 0.3083333 0.0000000 0.0000000
#> 
#> $pop1.pop3
#> # A folded site frequency spectrum of 2 populations
#> # Allele copies: pop1 = 2, pop3 = 2
#> # Loci used: 6 (dropped for missing data: 0)
#>     pop3
#> pop1         0   1         2
#>    0 0.7666667 1.8 0.3833333
#>    1 1.4666667 1.2 0.0000000
#>    2 0.3833333 0.0 0.0000000
#> 
#> $pop2.pop3
#> # A folded site frequency spectrum of 2 populations
#> # Allele copies: pop2 = 2, pop3 = 2
#> # Loci used: 6 (dropped for missing data: 0)
#>     pop3
#> pop2    0 1    2
#>    0 1.50 2 0.25
#>    1 1.00 1 0.00
#>    2 0.25 0 0.00
#> 
```
