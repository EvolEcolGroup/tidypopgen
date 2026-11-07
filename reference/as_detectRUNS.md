# Convert runs of homozygosity to the format of `detectRUNS`

Converts the runs found by
[`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md)
to the data frame used by the `detectRUNS` package, so that its
functions to plot and summarise runs (e.g.
[`detectRUNS::plot_Runs()`](https://rdrr.io/pkg/detectRUNS/man/plot_Runs.html),
[`detectRUNS::summaryRuns()`](https://rdrr.io/pkg/detectRUNS/man/summaryRuns.html))
can be used.

## Usage

``` r
as_detectRUNS(x)
```

## Arguments

- x:

  a tibble of runs, as returned by
  [`windows_indiv_roh()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_indiv_roh.md)

## Value

a data frame with columns `group` (the grouping variables, pasted with
"." if there are several, or the `id` if the runs were found on an
ungrouped
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)),
`id`, `chrom`, `nSNP`, `from`, `to` and `lengthBps`, in the order of `x`

## Examples

``` r
bed_file <- system.file("extdata", "roh", "roh_sim.bed",
  package = "tidypopgen"
)
sim_gt <- gen_tibble(bed_file, backingfile = tempfile(), quiet = TRUE)
sim_roh <- windows_indiv_roh(sim_gt %>% group_by(population))
head(as_detectRUNS(sim_roh))
#>      group          id chrom nSNP    from      to lengthBps
#> 1 moderate moderate_01     1  242 6221737 8461767   2240031
#> 2 moderate moderate_02     2  181 4990858 6269528   1278671
#> 3 moderate moderate_03     1  217 2355826 4534179   2178354
#> 4 moderate moderate_04     1  293 1244824 4160156   2915333
#> 5 moderate moderate_05     2  240 5845283 7391314   1546032
#> 6 moderate moderate_06     3  238 2703613 4112239   1408627
```
