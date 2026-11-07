# Block bootstrap of site frequency spectra

Compute bootstrap replicates of the site frequency spectra of
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_pop_sfs.md)
and
[`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/nwise_pop_sfs.md),
by resampling blocks of loci with replacement. Replicates are used to
estimate the uncertainty of statistics computed from the spectra, or of
the parameters of demographic models fitted with `dadi`, `moments` or
`fastsimcoal2` (see
[`write_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/write_sfs.md)
to write the replicates in their formats).

## Usage

``` r
bootstrap_pop_sfs(.x, ...)

# S3 method for class 'gen_tbl'
bootstrap_pop_sfs(
  .x,
  n_bootstrap = 100,
  block_size,
  size_unit = c("snp", "cM", "bp"),
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  ...
)

bootstrap_pairwise_pop_sfs(.x, ...)

# S3 method for class 'gen_tbl'
bootstrap_pairwise_pop_sfs(
  .x,
  n_bootstrap = 100,
  block_size,
  size_unit = c("snp", "cM", "bp"),
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
  n_cores = bigstatsr::nb_cores(),
  ...
)

bootstrap_nwise_pop_sfs(.x, ...)

# S3 method for class 'gen_tbl'
bootstrap_nwise_pop_sfs(
  .x,
  n_bootstrap = 100,
  block_size,
  size_unit = c("snp", "cM", "bp"),
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  fastest = c("last", "first"),
  max_cells = 1e+07,
  n_cores = bigstatsr::nb_cores(),
  ...
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md),
  grouped by population (or, for `bootstrap_pop_sfs()`, ungrouped to
  compute the spectrum of all individuals).

- ...:

  other arguments passed to specific methods, currently unused.

- n_bootstrap:

  the number of bootstrap replicates.

- block_size:

  the size of the blocks of loci that are resampled.

- size_unit:

  the unit of `block_size`: "snp" (the default), "cM" or "bp".

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

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- fastest:

  the order of the cells of each spectrum in the tidy output: "last"
  (the default) keeps the order of the populations, so that the rows are
  sorted by the number of allele copies in the first population, then in
  the second (the number of copies in the second population changes
  fastest, as in the files written by `dadi` and `easySFS`); "first"
  makes the number of copies in the first population change fastest.

- max_cells:

  the maximum number of cells of all the replicates (the number of cells
  of the joint spectrum times `n_bootstrap`); an error is raised for
  more cells, to avoid running out of memory.

## Value

If `type = "tidy"`, a tibble with the columns of the tidy output of
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_pop_sfs.md)
or
[`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/nwise_pop_sfs.md),
and a column `replicate` (after the group columns) with the number of
the replicate. If `type = "array"`, a list with one element per
replicate, each of which is what
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_pop_sfs.md)
or
[`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/nwise_pop_sfs.md)
return with `type = "array"` (lists of `sfs` objects, or an `sfs`
object). In both cases, the attribute `bootstrap_blocks` is a matrix
with the number of times each block (columns, numbered in the order in
which they appear in the loci table) was drawn in each replicate (rows),
and the tibble also has the attributes `n_proj` and `folded`.

## Details

The loci are divided into non-overlapping blocks that never span more
than one chromosome, defined as in the `jackknife_*` functions (e.g.
[`jackknife_pop_tajimas_d()`](https://evolecolgroup.github.io/tidypopgen/reference/jackknife_pop_tajimas_d.md)):
blocks of `block_size` SNPs, base pairs or centiMorgans (the latter
requires the `genetic_dist` column of the loci table). Blocks should be
larger than the extent of linkage disequilibrium, so that they are
approximately independent. Each replicate draws, with replacement, as
many blocks as there are (as the bootstrap of `dadi`,
`Misc.bootstraps_from_dd_chunks()`). With `block_size = 1` and
`size_unit = "snp"`, single SNPs are resampled (as the `--resample-vcf`
option of `easySFS`), which ignores linkage.

A spectrum is a sum over loci, so the spectrum of a replicate is the sum
of the spectra of the blocks it drew (counting each block as many times
as it was drawn). The genotypes are therefore read only once, to compute
the spectrum of each block. Missing data are dealt with by projection,
as in
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
with the same projections (`n_proj`) for all replicates: if `n_proj` is
`NULL`, they are chosen from all the loci with
[`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md).
The blocks are drawn with the random number generator of R, so use
[`set.seed()`](https://rdrr.io/r/base/Random.html) to make the
replicates reproducible.

## See also

[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs.md),
[`write_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/write_sfs.md),
[`jackknife_stats_generic()`](https://evolecolgroup.github.io/tidypopgen/reference/jackknife_stats_generic.md)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

set.seed(1)
# resample blocks of 2 SNPs
boot <- example_gt %>%
  bootstrap_pop_sfs(n_bootstrap = 10, block_size = 2, n_proj = 2)
boot
#> # A tibble: 90 × 5
#>    group replicate allele_count n_sites masked
#>    <chr>     <int>        <int>   <dbl> <lgl> 
#>  1 pop1          1            0    3.47 TRUE  
#>  2 pop1          1            1    2.53 FALSE 
#>  3 pop1          1            2    0    TRUE  
#>  4 pop1          2            0    3.6  TRUE  
#>  5 pop1          2            1    2.4  FALSE 
#>  6 pop1          2            2    0    TRUE  
#>  7 pop1          3            0    4    TRUE  
#>  8 pop1          3            1    2    FALSE 
#>  9 pop1          3            2    0    TRUE  
#> 10 pop1          4            0    3.2  TRUE  
#> # ℹ 80 more rows

# the 95% interval of the number of sites in each cell
boot %>%
  dplyr::group_by(group, allele_count) %>%
  dplyr::summarise(
    lower = stats::quantile(n_sites, 0.025),
    upper = stats::quantile(n_sites, 0.975)
  )
#> `summarise()` has regrouped the output.
#> ℹ Summaries were computed grouped by group and allele_count.
#> ℹ Output is grouped by group.
#> ℹ Use `summarise(.groups = "drop_last")` to silence this message.
#> ℹ Use `summarise(.by = c(group, allele_count))` for per-operation grouping
#>   (`?dplyr::dplyr_by`) instead.
#> # A tibble: 9 × 4
#> # Groups:   group [3]
#>   group allele_count lower upper
#>   <chr>        <int> <dbl> <dbl>
#> 1 pop1             0  2.96  3.91
#> 2 pop1             1  2.09  3.04
#> 3 pop1             2  0     0   
#> 4 pop2             0  3.5   4.5 
#> 5 pop2             1  1.5   2.5 
#> 6 pop2             2  0     0   
#> 7 pop3             0  1.5   4.5 
#> 8 pop3             1  1.5   4.5 
#> 9 pop3             2  0     0   

# the joint spectra of each pair, as arrays
boot_2d <- example_gt %>%
  bootstrap_pairwise_pop_sfs(
    n_bootstrap = 5, block_size = 2, n_proj = 2,
    type = "array"
  )
boot_2d[[1]]$pop1.pop2
#> # A folded site frequency spectrum of 2 populations
#> # Allele copies: pop1 = 2, pop2 = 2
#> # Loci used: 6 (dropped for missing data: 0)
#>     pop2
#> pop1         0 1         2
#>    0 0.9666667 1 0.4833333
#>    1 2.0666667 1 0.0000000
#>    2 0.4833333 0 0.0000000
```
