# Islands of runs of homozygosity

Finds the regions of the genome where many individuals (at least a
proportion `threshold` of them, overall or in each group) are in runs of
homozygosity: islands of runs of homozygosity (also called ROH
hotspots), which can point to regions under selection.

## Usage

``` r
gt_roh_islands(roh, loci, threshold = 0.5, by_group = TRUE)
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

- threshold:

  the smallest proportion of individuals in a run for a locus to be in
  an island (default 0.5)

- by_group:

  if `roh` was found on a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  whether to count the individuals of each group (the default) or of all
  groups

## Value

a tibble with one row per island, with the grouping variables (if
counting by group), `chromosome`, `start` and `end` (the positions of
its first and last loci), `length_bp` (`end - start + 1`), `n_loci`,
`first_locus` and `last_locus` (rows of the loci table), and `max_prop`
and `mean_prop` (the largest and mean proportions of individuals in a
run over its loci)

## Details

Islands are stretches of consecutive loci of a chromosome that are in
runs in at least `threshold` of the individuals (see
[`gt_roh_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_roh_loci.md)).

## Examples

``` r
bed_file <- system.file("extdata", "roh", "roh_sim.bed",
  package = "tidypopgen"
)
sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
# regions in runs in at least 20% of the individuals of a population
gt_roh_islands(sim_roh, show_loci(sim_gt), threshold = 0.2)
#> # A tibble: 8 × 10
#>   group    chromosome   start     end length_bp n_loci first_locus last_locus
#>   <chr>    <fct>        <int>   <int>     <int>  <int>       <int>      <int>
#> 1 inbred   1          1714516 8848406   7133891    724         165        888
#> 2 inbred   2          4700515 7489836   2789322    424        1507       1930
#> 3 inbred   3           754836 5045384   4290549    732        2131       2862
#> 4 moderate 1          2355826 4160156   1804331    179         227        405
#> 5 moderate 2          5845283 6269528    424246     67        1682       1748
#> 6 moderate 3          1495251 1802415    307165     59        2244       2302
#> 7 moderate 3          2703613 3009866    306254     63        2456       2518
#> 8 moderate 3          3763235 4112239    349005     53        2641       2693
#> # ℹ 2 more variables: max_prop <dbl>, mean_prop <dbl>
```
