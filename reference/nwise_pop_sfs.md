# Compute the joint site frequency spectrum of all populations

Compute the multidimensional (joint) site frequency spectrum (SFS) of
all the populations of a grouped `gen_tibble`, i.e. the number of loci
with each combination of numbers of copies of the derived (or minor)
allele in the populations (the "multiSFS" of `easySFS` and
`fastsimcoal2`). Missing genotypes are dealt with as in `easySFS`, by
projecting each population down to `n_proj` allele copies (see
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)
for details). Only the loci with at least `n_proj` called allele copies
in every population are used. A folded spectrum is folded on the total
number of copies across all populations (as in `dadi`).

## Usage

``` r
nwise_pop_sfs(.x, ...)

# S3 method for class 'gen_tbl'
nwise_pop_sfs(
  .x,
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
  max_cells = 1e+07,
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

  the order of the cells in the tidy output: "last" (the default) keeps
  the order of the populations, so that the rows are sorted by the
  number of allele copies in the first population, then in the second,
  and so on (the number of copies in the last population changes
  fastest, as in the files written by `dadi` and `easySFS`); "first"
  makes the number of copies in the first population change fastest.

- max_cells:

  the maximum number of cells of the spectrum; an error is raised for
  larger spectra, to avoid running out of memory.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- block_size:

  maximum number of loci read at once.

## Value

If `type = "tidy"`, a tibble with one row per cell of the joint
spectrum, with one column per population (named with the group labels)
giving its number of copies of the derived (or minor) allele (see
`allele_count` in
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)),
`n_sites` (the projected number of sites) and `masked` (whether the cell
is masked, see
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)).
The projection of each population is stored in the attribute `n_proj`.
If `type = "array"`, an `sfs` object (see
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md)),
an array with one dimension per population.

## Details

Populations are in the order of the groups (see
[`dplyr::group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html)):
to change their order (which matters for `dadi` and `fastsimcoal2`),
change the order of the levels of the grouping variable (e.g. with
[`factor()`](https://rdrr.io/r/base/factor.html)); to use a subset of
the populations,
[filter()](https://dplyr.tidyverse.org/reference/filter.html) the
`gen_tibble` (empty groups are dropped). The spectrum has
`prod(n_proj + 1)` cells, which grows quickly with the number of
populations.

## See also

[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_pop_sfs.md),
[`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

example_gt %>% nwise_pop_sfs(n_proj = 2)
#> # A tibble: 27 × 5
#>     pop1  pop2  pop3 n_sites masked
#>    <int> <int> <int>   <dbl> <lgl> 
#>  1     0     0     0  0.617  TRUE  
#>  2     0     0     1  0.967  FALSE 
#>  3     0     0     2  0.200  FALSE 
#>  4     0     1     0  0.117  FALSE 
#>  5     0     1     1  0.533  FALSE 
#>  6     0     1     2  0.142  FALSE 
#>  7     0     2     0  0.0333 FALSE 
#>  8     0     2     1  0.15   FALSE 
#>  9     0     2     2  0      TRUE  
#> 10     1     0     0  0.6    FALSE 
#> # ℹ 17 more rows

# the first population changing fastest
example_gt %>% nwise_pop_sfs(n_proj = 2, fastest = "first")
#> # A tibble: 27 × 5
#>     pop1  pop2  pop3 n_sites masked
#>    <int> <int> <int>   <dbl> <lgl> 
#>  1     0     0     0  0.617  TRUE  
#>  2     1     0     0  0.6    FALSE 
#>  3     2     0     0  0.283  FALSE 
#>  4     0     1     0  0.117  FALSE 
#>  5     1     1     0  0.6    FALSE 
#>  6     2     1     0  0.142  FALSE 
#>  7     0     2     0  0.0333 FALSE 
#>  8     1     2     0  0.133  FALSE 
#>  9     2     2     0  0      TRUE  
#> 10     0     0     1  0.967  FALSE 
#> # ℹ 17 more rows
```
