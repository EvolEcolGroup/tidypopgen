# Tidy a site frequency spectrum

Turn an `sfs` object (as returned by
[`pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pop_sfs.md),
[`pairwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_sfs.md)
and
[`nwise_pop_sfs()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/nwise_pop_sfs.md)
with `type = "array"`) into a tibble with one row per cell of the
spectrum.

## Usage

``` r
# S3 method for class 'sfs'
tidy(x, fastest = c("last", "first"), ...)
```

## Arguments

- x:

  an `sfs` object

- fastest:

  which population's number of allele copies changes fastest from one
  row to the next: "last" (the default, which keeps the order of the
  populations: rows are sorted by the number of copies in the first
  population, then in the second, and so on, as in the files written by
  dadi and easySFS) or "first".

- ...:

  currently unused.

## Value

a tibble with one column per population (named with the population
labels) giving the number of allele copies (for an unfolded spectrum,
copies of the derived allele: the alternate allele, unless the loci were
polarised with `ancestral`, in which case it is the reference allele at
loci whose ancestral allele is the alternate one; for a folded spectrum,
copies of the minor allele across all the populations of the spectrum),
`n_sites`, the number of sites in that cell (a fractional number, as
each locus is spread over the cells by the projection), and `masked`,
TRUE for the cells that are masked (the corners, i.e. monomorphic sites,
and for a folded spectrum the cells that were folded out). The
projection of each population is stored in the attribute `n_proj`, and
whether the spectrum is folded in the attribute `folded`.

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

sfs <- example_gt %>% pop_sfs(n_proj = 4, type = "array")
tidy(sfs[[1]])
#> # A tibble: 5 × 3
#>    pop1 n_sites masked
#>   <int>   <dbl> <lgl> 
#> 1     0    1.47 TRUE  
#> 2     1    2.13 FALSE 
#> 3     2    2.4  FALSE 
#> 4     3    0    TRUE  
#> 5     4    0    TRUE  
```
