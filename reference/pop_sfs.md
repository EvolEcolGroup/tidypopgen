# Compute the site frequency spectrum of each population

Compute the one-dimensional site frequency spectrum (SFS) of each
population, i.e. the number of loci with 0, 1, 2, ... copies of the
derived (or minor) allele. Missing genotypes are dealt with as in
`easySFS` (which relies on `dadi`), so that the spectra can be used for
demographic inference with `dadi`, `moments` or `fastsimcoal2`. See
[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_pop_sfs.md)
for the joint spectra of pairs of populations, and
[`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/reference/nwise_pop_sfs.md)
for the joint spectrum of all populations.

## Usage

``` r
pop_sfs(.x, ...)

# S3 method for class 'tbl_df'
pop_sfs(
  .x,
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
pop_sfs(
  .x,
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

# S3 method for class 'grouped_gen_tbl'
pop_sfs(
  .x,
  n_proj = NULL,
  folded = TRUE,
  ancestral = NULL,
  type = c("tidy", "array"),
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  (grouped by population, or ungrouped to compute a single spectrum for
  all individuals), or a vector of class `vctrs_bigSNP` (usually the
  `genotypes` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)).

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

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- block_size:

  maximum number of loci read at once.

## Value

If `type = "tidy"`, a tibble with one row per cell of the spectrum of
each population, with the columns:

- `group` (or one column per grouping variable, if the data are grouped
  by more than one variable), absent for ungrouped data;

- `allele_count`: the number of copies of the derived allele (for an
  unfolded spectrum: the alternate allele, unless the loci are polarised
  with `ancestral`) or of the minor allele (for a folded spectrum);

- `n_sites`: the (projected) number of sites;

- `masked`: whether the cell is masked (see details).

The projection of each population is stored in the attribute `n_proj`.
The tibble has class `pop_sfs_tbl`, and can be plotted with
[autoplot()](https://evolecolgroup.github.io/tidypopgen/reference/autoplot_sfs.md).
If `type = "array"`, a list of `sfs` objects (one per population, named
by the group labels), or a single `sfs` object for ungrouped data. An
`sfs` object is an array of the numbers of sites, with attributes giving
its `mask`, the population labels (`pop_ids`), `n_proj`, whether it is
`folded`, and the number of loci used (`n_loci`), dropped for lack of
called allele copies (`n_loci_dropped`), and dropped because their
ancestral allele is unknown (`n_loci_unpolarised`, only with
`ancestral`). Use
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) to turn it
into a tibble.

## Details

**Missing data.** Because of missing genotypes, the number of allele
copies called at a locus varies from locus to locus. Each population is
therefore projected down to a fixed number of allele copies, `n_proj`
(Marth et al. 2004; Gutenkunst et al. 2009): a locus with \\n\\ called
allele copies, \\k\\ of which are derived, contributes to cell \\j\\ of
the spectrum the hypergeometric probability of drawing \\j\\ derived
copies when sampling `n_proj` of its \\n\\ copies without replacement.
Loci with fewer than `n_proj` called copies are dropped. Each locus is
thus spread over several cells, and the spectrum contains fractional
numbers of sites. This is the approach of `easySFS` (with the `-a`
option, i.e. using all loci).

**Choosing the projection.** Projecting to more allele copies gives more
resolution (more segregating sites, especially rare variants), but drops
more loci with missing data. Following the advice of the `dadi` manual
and `easySFS`, use
[`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md)
to compute the number of segregating sites for each possible projection,
and choose a value that keeps most segregating sites. If `n_proj` is not
given, each population is projected to the number of allele copies that
maximises its number of segregating sites (the `best` projection in
[`pop_sfs_preview()`](https://evolecolgroup.github.io/tidypopgen/reference/pop_sfs_preview.md)).

**Folding and polarisation.** Unless `folded = FALSE`, the spectrum is
folded (i.e. it counts copies of the minor allele), as the ancestral
allele is not usually known. An unfolded spectrum counts copies of the
derived allele. By default (`ancestral = NULL`), the reference allele is
taken to be ancestral, as in `easySFS`, so that the spectrum counts
copies of the alternate allele: this is only meaningful if the alleles
of the `gen_tibble` are polarised (note that some operations, such as
merging datasets with [`rbind()`](https://rdrr.io/r/base/cbind.html),
can swap the reference and alternate alleles). Otherwise, `ancestral`
gives the ancestral allele of each locus (e.g. from an outgroup, or the
`AA` field of a vcf), either as a column of the loci table or as a
vector: loci whose ancestral allele is the alternate allele are
polarised (copies of the reference allele are counted), and loci whose
ancestral allele is missing (`NA`) or is neither of their two alleles
are dropped, as in `dadi`. Alleles are compared ignoring case.
Polarisation does not change a folded spectrum, so `ancestral` is
ignored (with a warning) when `folded = TRUE`.

**Masked cells.** As in `dadi`, the corners of the spectrum (sites
monomorphic in the projected sample) are masked, and so are the cells
that were folded out in a folded spectrum. They are reported in the
`masked` column (or the `mask` attribute), but keep their values.

Each spectrum only uses the loci that have enough called allele copies
in its own populations. Diploid, polyploid, mixed-ploidy and
pseudohaploid data are supported: each individual contributes as many
allele copies as its ploidy (one for pseudohaploids). The spectrum is
computed from the genotypes currently in use: if the `gen_tibble` uses
imputed genotypes (see
[`gt_set_imputed()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_set_imputed.md)),
these are treated as called genotypes, so you will usually want to use
`gt_set_imputed(x, FALSE)` first.

Note that `easySFS` samples a single SNP per `CHROM` (i.e. per RAD
locus) by default; to do the same, thin the loci (e.g. with
[`select_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/select_loci.md))
before computing the spectra.

## References

Marth GT, Czabarka E, Murvai J, Sherry ST (2004) The allele frequency
spectrum in genome-wide human variation data reveals signals of
differential demographic history in three large world populations.
Genetics 166: 351-372.

Gutenkunst RN, Hernandez RD, Williamson SH, Bustamante CD (2009)
Inferring the joint demographic history of multiple populations from
multidimensional SNP frequency data. PLoS Genetics 5: e1000695.

Overcast I. easySFS. <https://github.com/isaacovercast/easySFS>

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Preview the number of segregating sites for each projection
example_gt %>% pop_sfs_preview()
#> # A tibble: 11 × 5
#>    group n_proj n_seg_sites n_loci best 
#>    <chr>  <int>       <dbl>  <dbl> <lgl>
#>  1 pop1       2        2.67      6 FALSE
#>  2 pop1       3        4         6 FALSE
#>  3 pop1       4        4.53      6 TRUE 
#>  4 pop1       5        3.83      5 FALSE
#>  5 pop1       6        4         5 FALSE
#>  6 pop2       2        2         6 FALSE
#>  7 pop2       3        3         5 FALSE
#>  8 pop2       4        4         5 TRUE 
#>  9 pop3       2        3         6 TRUE 
#> 10 pop3       3        1.5       4 FALSE
#> 11 pop3       4        2         4 FALSE

# Compute the folded spectra, projecting each population to 2 copies
example_gt %>% pop_sfs(n_proj = 2)
#> # A tibble: 9 × 4
#>   group allele_count n_sites masked
#>   <chr>        <int>   <dbl> <lgl> 
#> 1 pop1             0    3.33 TRUE  
#> 2 pop1             1    2.67 FALSE 
#> 3 pop1             2    0    TRUE  
#> 4 pop2             0    4    TRUE  
#> 5 pop2             1    2    FALSE 
#> 6 pop2             2    0    TRUE  
#> 7 pop3             0    3    TRUE  
#> 8 pop3             1    3    FALSE 
#> 9 pop3             2    0    TRUE  

# Project each population to a different number of copies
example_gt %>% pop_sfs(n_proj = c(pop1 = 4, pop2 = 2, pop3 = 2))
#> # A tibble: 11 × 4
#>    group allele_count n_sites masked
#>    <chr>        <int>   <dbl> <lgl> 
#>  1 pop1             0    1.47 TRUE  
#>  2 pop1             1    2.13 FALSE 
#>  3 pop1             2    2.4  FALSE 
#>  4 pop1             3    0    TRUE  
#>  5 pop1             4    0    TRUE  
#>  6 pop2             0    4    TRUE  
#>  7 pop2             1    2    FALSE 
#>  8 pop2             2    0    TRUE  
#>  9 pop3             0    3    TRUE  
#> 10 pop3             1    3    FALSE 
#> 11 pop3             2    0    TRUE  

# Unfolded spectra, as arrays
example_gt %>% pop_sfs(n_proj = 2, folded = FALSE, type = "array")
#> $pop1
#> # An unfolded site frequency spectrum of 1 population
#> # Allele copies: pop1 = 2
#> # Loci used: 6 (dropped for missing data: 0)
#> pop1
#>        0        1        2 
#> 2.166667 2.666667 1.166667 
#> 
#> $pop2
#> # An unfolded site frequency spectrum of 1 population
#> # Allele copies: pop2 = 2
#> # Loci used: 6 (dropped for missing data: 0)
#> pop2
#>   0   1   2 
#> 3.5 2.0 0.5 
#> 
#> $pop3
#> # An unfolded site frequency spectrum of 1 population
#> # Allele copies: pop3 = 2
#> # Loci used: 6 (dropped for missing data: 0)
#> pop3
#> 0 1 2 
#> 3 3 0 
#> 

# Unfolded spectra, polarised with the ancestral allele of each locus (the
# locus with an unknown ancestral allele is dropped)
show_loci(example_gt)
#> # A tibble: 6 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs1   chr1              3            0 A          T         
#> 2         2 rs2   chr1              5            0 T          C         
#> 3         3 rs3   chr1             65            0 C          NA        
#> 4         4 rs4   chr1            343            0 G          C         
#> 5         5 rs5   chr2             23            0 C          G         
#> 6         6 rs6   chr2            456            0 T          A         
example_gt %>%
  pop_sfs(
    n_proj = 2, folded = FALSE,
    ancestral = c("T", "C", "C", "G", NA, "T")
  )
#> # A tibble: 9 × 4
#>   group allele_count n_sites masked
#>   <chr>        <int>   <dbl> <lgl> 
#> 1 pop1             0    2.63 TRUE  
#> 2 pop1             1    2.07 FALSE 
#> 3 pop1             2    0.3  TRUE  
#> 4 pop2             0    2.5  TRUE  
#> 5 pop2             1    1.5  FALSE 
#> 6 pop2             2    1    TRUE  
#> 7 pop3             0    1.5  TRUE  
#> 8 pop3             1    2    FALSE 
#> 9 pop3             2    1.5  TRUE  
```
