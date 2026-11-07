# Inbreeding from runs of homozygosity (F_ROH)

Computes F_ROH, the proportion of the genome of each individual that is
in runs of homozygosity, from the runs found by
[`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md).
F_ROH can be computed from the runs of at least a given length (e.g. 1,
2, 4, 8 and 16 Mb): long runs come from recent common ancestors, while
short runs also reflect older inbreeding.

## Usage

``` r
gt_roh_froh(
  roh,
  min_length_mb = 0,
  by = c("genome", "chromosome"),
  genome_length = NULL
)
```

## Arguments

- roh:

  an `indiv_roh` object, as returned by
  [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md)

- min_length_mb:

  the smallest length of the runs to use, in Mb (default 0, i.e. all
  runs). With several values, F_ROH is computed for each.

- by:

  whether to compute F_ROH over the whole genome (`genome`, the
  default), or for each chromosome (`chromosome`)

- genome_length:

  the length of the genome, in bp (only used with `by = "genome"`).
  Default is the sum of the lengths of the chromosomes (see Details).

## Value

a tibble with one row per individual (including those with no runs),
with the grouping variables of `roh` (if any), `id`, the column
`min_length_mb` if there are several values, `chromosome` if
`by = "chromosome"`, `length_bp` (the total length of the runs) and
`froh`.

## Details

The length of the genome is, by default, the sum of the lengths of the
chromosomes, each being its largest position (as in `detectRUNS`); it
can be given instead with `genome_length` (e.g. the length of the
autosomes of the reference genome, or the length covered by the loci).

## See also

[`autoplot.indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/autoplot.indiv_roh.md)
to plot F_ROH.

## Examples

``` r
bed_file <- system.file("extdata", "roh", "roh_sim.bed",
  package = "tidypopgen"
)
sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
# F_ROH of each individual
gt_roh_froh(sim_roh)
#> # A tibble: 30 × 4
#>    group   id         length_bp  froh
#>    <chr>   <chr>          <dbl> <dbl>
#>  1 outbred outbred_01         0     0
#>  2 outbred outbred_02         0     0
#>  3 outbred outbred_03         0     0
#>  4 outbred outbred_04         0     0
#>  5 outbred outbred_05         0     0
#>  6 outbred outbred_06         0     0
#>  7 outbred outbred_07         0     0
#>  8 outbred outbred_08         0     0
#>  9 outbred outbred_09         0     0
#> 10 outbred outbred_10         0     0
#> # ℹ 20 more rows
# from the runs of at least 1, 2 and 4 Mb
gt_roh_froh(sim_roh, min_length_mb = c(1, 2, 4))
#> # A tibble: 90 × 5
#>    group   id         min_length_mb length_bp  froh
#>    <chr>   <chr>              <dbl>     <dbl> <dbl>
#>  1 outbred outbred_01             1         0     0
#>  2 outbred outbred_01             2         0     0
#>  3 outbred outbred_01             4         0     0
#>  4 outbred outbred_02             1         0     0
#>  5 outbred outbred_02             2         0     0
#>  6 outbred outbred_02             4         0     0
#>  7 outbred outbred_03             1         0     0
#>  8 outbred outbred_03             2         0     0
#>  9 outbred outbred_03             4         0     0
#> 10 outbred outbred_04             1         0     0
#> # ℹ 80 more rows
# on each chromosome
gt_roh_froh(sim_roh, by = "chromosome")
#> # A tibble: 90 × 5
#>    group   id         chromosome length_bp  froh
#>    <chr>   <chr>      <fct>          <dbl> <dbl>
#>  1 outbred outbred_01 1                  0     0
#>  2 outbred outbred_01 2                  0     0
#>  3 outbred outbred_01 3                  0     0
#>  4 outbred outbred_02 1                  0     0
#>  5 outbred outbred_02 2                  0     0
#>  6 outbred outbred_02 3                  0     0
#>  7 outbred outbred_03 1                  0     0
#>  8 outbred outbred_03 2                  0     0
#>  9 outbred outbred_03 3                  0     0
#> 10 outbred outbred_04 1                  0     0
#> # ℹ 80 more rows
```
