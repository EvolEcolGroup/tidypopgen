# The loci in runs of homozygosity

Counts, for each locus, the individuals that are in a run of
homozygosity (as found by
[`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md)),
overall or in each group. Loci that are in runs in most individuals of a
population (see
[`gt_roh_islands()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_islands.md))
can point to regions under selection. The counts are those of the
`.hom.summary` file of PLINK `--homozyg`.

## Usage

``` r
gt_roh_loci(roh, loci, by_group = TRUE)
```

## Arguments

- roh:

  an `indiv_roh` object, as returned by
  [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/windows_indiv_roh.md)

- loci:

  the loci table of the
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  used to find the runs (i.e. `show_loci(x)`, as the runs refer to its
  rows). Integer chromosomes (`chromosomes_as_int = TRUE`) are replaced
  by the names of the chromosomes of the runs, so the output is the same
  as without `chromosomes_as_int`; integer names of loci
  (`names_as_int = TRUE`) are kept (see
  [`gt_restore_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_lookup_loci.md)
  to restore them first)

- by_group:

  if `roh` was found on a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  whether to count the individuals of each group (the default) or of all
  groups

## Value

a tibble with one row per locus (and per group), with the grouping
variables (if counting by group), `name`, `chromosome`, `position`,
`n_roh` (the number of individuals in a run), `n_indiv` (the number of
individuals) and `prop_roh` (`n_roh / n_indiv`)

## See also

[`gt_roh_islands()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_islands.md)
to find the regions in runs in many individuals, and
[`autoplot.indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/autoplot.indiv_roh.md)
(`type = "coverage"`) to plot them.

## Examples

``` r
bed_file <- system.file("extdata", "roh", "roh_sim.bed",
  package = "tidypopgen"
)
sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
roh_loci <- gt_roh_loci(sim_roh, show_loci(sim_gt))
# the loci most often in runs in each population
roh_loci %>%
  group_by(group) %>%
  slice_max(prop_roh, n = 1, with_ties = FALSE)
#> # A tibble: 3 × 7
#> # Groups:   group [3]
#>   group    name    chromosome position n_roh n_indiv prop_roh
#>   <chr>    <chr>   <fct>         <int> <int>   <int>    <dbl>
#> 1 inbred   snp2551 3           3230132     7      10      0.7
#> 2 moderate snp227  1           2355826     2      10      0.2
#> 3 outbred  snp1    1             24153     0      10      0  
```
